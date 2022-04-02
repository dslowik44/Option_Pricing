#include "option_pricing.hpp"
// #include "fexcpt.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <iostream>
#include <random>

using std::exp;
using std::sqrt;

namespace OptionPricing {

    Asset::Asset(double S_o, double q, double sigma)
        : _S_o{ S_o }, _q{ q }, _sigma{ sigma } { }


    Terms::Terms(Style s, Type t, double T, double K, AvgType avg)
        // Set terms of Euro, Amer, or AsianPrc Call or Put, each of which define their 
        // _payoff(S) in terms of the strike K in the usual vanilla way.
        try : _style{ s }, _type{ t }, _T{ T }, _K{ K }, _avg_type{ avg } {
        if (t == Type::Other || s == Style::AsianStrike || s == Style::AsianExotic
            || avg == AvgType::Other || _K < 0 || _T <= 0) {
            throw std::invalid_argument(
                "Must provide Euro, Amer or AsianPrc Call or Put with valid K and T.");
        }
        if (s == Style::AsianPrc && avg != AvgType::Arith && avg != AvgType::Geom)
            throw std::invalid_argument("AsianPrc must do an Arith or Geom path price average.\n \
                                         -Use AsianExotic to specify Other path price avgeraging.");
        if (_type == Type::Call)
            _payoff = [K](double S) { return std::max(S - K, 0.0); };
        else
            _payoff = [K](double S) { return std::max(K - S, 0.0); };
        if (_style == Style::AsianPrc) {
            if (_avg_type == AvgType::Arith) {
                _path_accum_fn = arith_path_accum_fn;
                _init_accum = 0.0;
            } else {
                _path_accum_fn = geom_path_accum_fn;
                _init_accum = 1.0;
            }
        } else  // Euro and Amer have no path averaging:
            _avg_type = AvgType::None;
    }
    catch (std::exception& e) {
        std::cout << "In Terms's ctor1: " << e.what() << std::endl;
        throw;
    }

    Terms::Terms(Style s, Type t, double T, AvgType avg)
        // Set terms of AsianStrike. Like AsianPrc, must be a Call or Put and uses only Arith or 
        // Geom path price averaging, -Use AsianExotic to specify other path averaging.
        try : _style{ s }, _type{ t }, _T{ T }, _K{ -1 }, _avg_type{ avg } {
        if (s != Style::AsianStrike || t == Type::Other ||
            (avg != AvgType::Arith && avg != AvgType::Geom) || _T <= 0) {
            throw std::invalid_argument("Must provide AsianStrike Call or Put with valid T\n \
             and can only do Arith or Geom path price average:\n \
              -Use AsianExotic to specify Other path price averaging.");
        }
        if (_type == Type::Call)
            _path_payoff = [](double K_avg, double S) { return std::max(S - K_avg, 0.0); };
        else
            _path_payoff = [](double K_avg, double S) { return std::max(K_avg - S, 0.0); };

        if (_avg_type == AvgType::Arith) {
            _path_accum_fn = arith_path_accum_fn;
            _init_accum = 0.0;
        } else {
            _path_accum_fn = geom_path_accum_fn;
            _init_accum = 1.0;
        }
    }
    catch (std::exception& e) {
        std::cout << "In Terms's ctor2: " << e.what() << std::endl;
        throw;
    }

    Terms::Terms(Style s, Type t, double T, std::function<double(double)> payoff)
        // Set terms of Euro or Amer option that is not Call or Put (i.e. is Other)
        try : _style{ s }, _type{ t }, _T{ T }, _K{ -1 }, _avg_type{ AvgType::None } {
        if ((s != Style::Amer && s != Style::Euro) || t != Type::Other || _T <= 0) {
            throw std::invalid_argument("Must provide Euro or Amer with valid T\n \
             and provide payoff function of this exotic(Type::Other) Terms object.");
        }
        _payoff = payoff;
    }
    catch (std::exception& e) {
        std::cout << "In Terms's ctor3: " << e.what() << std::endl;
        throw;
    }

    Terms::Terms(Style s, Type t, double T, std::function<double(double, double)> path_payoff,
        std::function<void(int, double, double&)> path_accum_fn, double init_accum)
        // Set terms of AsianExotic option. 
        // Both path_payoff and path_accum_fn (along with init_accum) need to be specified.
        try : _style{ s }, _type{ t }, _T{ T }, _K{ -1 }, _avg_type{ AvgType::Other } {
        if (s != Style::AsianExotic || t != Type::Other || _T <= 0) {
            throw std::invalid_argument("Must provide AsianExotic with valid T\n \
             path_payoff function of this exotic(Type::Other) Terms object.");
        }
        _path_payoff = path_payoff;
        _path_accum_fn = path_accum_fn;
        _init_accum = init_accum;
    }
    catch (std::exception& e) {
        std::cout << "In Terms's ctor4: " << e.what() << std::endl;
        throw;
    }


    std::function<void(int, double, double&)>
        Terms::arith_path_accum_fn { [](int j, double s, double& cum) {
        cum += (s - cum) / j;
    } };
    std::function<void(int, double, double&)>
        Terms::geom_path_accum_fn = [](int j, double s, double& cum) {
        cum *= pow(s / cum, 1.0 / j);
    };

    Option::Option(const Terms& t, const Asset& a)
        : _terms{ t }, _asset{ a }, _prc{ -1 }, _greeks{} {  }


    double Option::btree_prc(int N_steps, double r, bool get_greeks) {
        Greeks* g_ptr{ get_greeks ? &_greeks : nullptr };
        _prc = btree_prc_terms_assets(_terms, _asset, N_steps, r, g_ptr);
        return _prc;
    }

    double btree_prc_opt(Option& opt, int N_steps, double r, bool greeks) {
        return opt.btree_prc(N_steps, r, greeks);
    }

    double btree_prc_terms_assets(const Terms& t, const Asset& a, int N_steps, double r,
        Option::Greeks* g_ptr) {
        if (t.style() != Terms::Style::Amer && t.style() != Terms::Style::Euro) {
            throw std::invalid_argument("In btree_prc_raw, Bad Terms::Style.");
        }
        return btree_prc_raw(N_steps, r, t.style(), t.T(), t.payoff_fn(),
            a.S_o(), a.q(), a.sigma(), g_ptr);
    }

    double btree_prc_raw(int N, double r,
        Terms::Style style, double T, std::function<double(double)> payoff,
        double S_o, double q, double sigma, Option::Greeks* g_ptr) {

        if (style != Terms::Style::Amer && style != Terms::Style::Euro) {
            throw std::invalid_argument("In btree_prc_raw, Bad Terms::Style.");
        }

        double dt = T / N;
        double u = exp(sigma * sqrt(dt));
        double u2 = u * u;
        double d = 1. / u;
        double a = exp((r - q) * dt);
        double p = (a - d) / (u - d);
        double df = exp(-r * dt);

        struct Node {
            double s{ 0 };
            double f{ 0 };
        };
        int len_nodes = N + 1;
        Node* nodes = new Node[len_nodes]{};

        // Set s and f for final time step:
        nodes[0].s = S_o * std::pow(u, -N);
        nodes[0].f = payoff(nodes[0].s);
        for (int j = 1; j < N + 1; ++j) {
            if (j == N / 2) {  // Take advantage to tie node's s value to known value at middles:
                if (N % 2)     // Odd N, we're one step below middle:
                    nodes[j].s = S_o / u;
                else           // Even N, we're at the middle, so s == S_o:
                    nodes[j].s = S_o;
            } else {
                nodes[j].s = nodes[j - 1].s * u2;
            }
            nodes[j].f = payoff(nodes[j].s);
        }

        std::function<double(int)> step_back_Payoff;
        if (style == Terms::Style::Euro)
            step_back_Payoff = [p, df, nodes](int idx) {
            return df * (p * nodes[idx + 1].f + (1 - p) * nodes[idx].f); };
        else
            step_back_Payoff = [p, df, payoff, nodes](int idx) {
            return std::max(payoff(nodes[idx].s), df * (p * nodes[idx + 1].f + (1 - p) * nodes[idx].f)); };

        // Backstep from time i=N-1 to i=0, filling in nodes j=0 to j=i for each:
        double s1, s2, s3, s4, s5;  // Needed to save some values at i=1, 2 for greeks calculations.
        double f1, f2, f3, f4, f5;
        for (int i = N - 1; i >= 0; --i) {
            for (int j = 0; j <= i; ++j) {
                if (j == i / 2) {  // As above, tie to known values.
                    if (i % 2)     // Potential speed up: break j-loop into before and after middle..
                        nodes[j].s = S_o / u;
                    else
                        nodes[j].s = S_o;
                } else {
                    nodes[j].s = nodes[j].s * u;
                }
                nodes[j].f = step_back_Payoff(j);  // flag_Payoffe("Underflow for fine grid: Ok", false);
// see https://stackoverflow.com/questions/46611633/can-the-floating-point-status-flag-fe-underflow-set-when-the-result-is-not-sub-n
            }
            if (i <= 2) {  // Save some values to calculate delta, gamma, theta:
                if (i == 1) {
                    s1 = nodes[0].s;
                    s2 = nodes[1].s;
                    f1 = nodes[0].f;
                    f2 = nodes[1].f;
                }
                if (i == 2) {
                    s3 = nodes[0].s;
                    s4 = nodes[1].s;
                    s5 = nodes[2].s;
                    f3 = nodes[0].f;
                    f4 = nodes[1].f;
                    f5 = nodes[2].f;
                }
            }
        }

        double f0 = nodes[0].f;

        if (g_ptr) {
            g_ptr->delta = (f2 - f1) / (s2 - s1);
            g_ptr->gamma = 2 * ((f5 - f4) / (s5 - s4) - (f4 - f3) / (s4 - s3)) / (s5 - s3);
            g_ptr->theta = (f4 - f0) / (2 * dt);

            double incr = .0001;
            g_ptr->rho = (btree_prc_raw(N, r + incr, style, T, payoff, S_o, q, sigma, nullptr) - f0) / incr;
            g_ptr->vega = (btree_prc_raw(N, r, style, T, payoff, S_o, q, sigma + incr, nullptr) - f0) / incr;
            //g_ptr->delta2 = (btree_prc_raw(N, r, style, T, F, S_o + incr, q, sigma, nullptr) - f0) / incr;
            //g_ptr->theta2 = (btree_prc_raw(N, r, style, T - incr, F, S_o, q, sigma, nullptr) - f0) / incr;
        }

        delete[] nodes;
        return f0;
    }


    double Option::MC_prc(int N_sims, int N_steps, double r, bool get_greeks, bool antithetic, double& sd_prc) {
        Greeks* g_ptr{ get_greeks ? &_greeks : nullptr };
        _prc = MC_prc_terms_assets(_terms, _asset, N_sims, N_steps, r, g_ptr, antithetic, sd_prc);
        return _prc;
    }

    double MC_prc_opt(Option& opt, int N_sims, int N_steps, double r, bool greeks, bool antithetic, double& sd_prc) {
        return opt.MC_prc(N_sims, N_steps, r, greeks, antithetic, sd_prc);
    }

    double MC_prc_terms_assets(const Terms& t, const Asset& a, int N_sims, int N_steps, double r,
        Option::Greeks* g_ptr, bool antithetic, double& sd_prc) {
        if (t.style() == Terms::Style::Amer) {
            throw std::invalid_argument("In MC_prc_terms_assets, Can't price Amer option via MC.");
        }
        return MC_prc_raw(N_sims, N_steps, r, t.style(), t.T(), t.payoff_fn(), t.path_payoff(),
            t.path_accum(), t.init_accum(), a.S_o(), a.q(), a.sigma(), g_ptr, antithetic, sd_prc);
    }

    double MC_prc_raw(int N_sims, int N_steps, double r, Terms::Style style, double T, std::function<double(double)> payoff,
        std::function<double(double, double)> path_payoff, std::function<void(int, double, double&)> path_accum_fn, double init_accum,
        double S_o, double q, double sigma, Option::Greeks* g_ptr, bool antithetic, double& sd_mean_PV) {

        if (style == Terms::Style::Amer) {
            throw std::invalid_argument("In MC_prc_raw, no Amer terms allowed. \n");
        }

        double g = r - q - sigma * sigma / 2;  // risk neutral growth rate of lnS.
        double df = exp(-r * T);

        // long unsigned int seed = 1;    // Replace next 2 lines with these 2 for reproducing random 
        // static std::mt19937 gen{seed}; // bit streams for testing.
        std::random_device rd;
        static std::mt19937 gen{ rd() };
        std::normal_distribution<> norm{ };

        double cumm_P{ 0.0 };
        double ssq{ 0.0 };
        if (style == Terms::Style::Euro) {
            // Can by-pass path stepping and just simulate final risk-neutral asset price at time T:
            double sd = sigma * sqrt(T);

            if (!g_ptr) {  // nullptr, skip greeks:
                if (!antithetic) {
                    for (int i = 0; i < N_sims; ++i) {
                        double eps = norm(gen);
                        double P = payoff(S_o * exp(g * T + eps * sd));
                        cumm_P += P;
                        ssq += P * P;
                    }
                } else {
                    for (int i = 0; i < N_sims; ++i) {
                        double eps = norm(gen);
                        double P = (payoff(S_o * exp(g * T + eps * sd)) +
                            payoff(S_o * exp(g * T - eps * sd))) / 2;
                        cumm_P += P;
                        ssq += P * P;
                    }
                }
            } else {      // valid g_ptr, get greeks:
                double incr = 0.01;
                double S_o_up{ (1. + incr) * S_o }, S_o_dn{ (1. - incr) * S_o }, g_rho_up{ g + incr * r };
                double cumm_delta_up{ 0 }, cumm_delta_dn{ 0 }, cumm_rho_up{ 0 }, cumm_theta_up{ 0 }, cumm_vega_up{ 0 };
                double sigma_up{ (1 + incr) * sigma };
                if (!antithetic) {
                    for (int i = 0; i < N_sims; ++i) {
                        double eps = norm(gen);
                        double P = payoff(S_o * exp(g * T + eps * sd));
                        cumm_P += P;
                        ssq += P * P;
                        cumm_delta_up += payoff(S_o_up * exp(g * T + eps * sd));
                        cumm_delta_dn += payoff(S_o_dn * exp(g * T + eps * sd));
                        cumm_rho_up += payoff(S_o * exp(g_rho_up * T + eps * sd));
                        cumm_theta_up += payoff(S_o * exp(g * (1 - incr) * T + eps * sigma * sqrt((1 - incr) * T)));
                        cumm_vega_up += payoff(S_o * exp((r - q - sigma_up * sigma_up / 2) * T + eps * sigma_up * sqrt(T)));
                    }
                } else {
                    for (int i = 0; i < N_sims; ++i) {
                        double eps = norm(gen);
                        double P = (payoff(S_o * exp(g * T + eps * sd)) +
                            payoff(S_o * exp(g * T - eps * sd))) / 2;
                        cumm_P += P;
                        ssq += P * P;
                        cumm_delta_up += (payoff(S_o_up * exp(g * T + eps * sd)) +
                            payoff(S_o_up * exp(g * T - eps * sd))) / 2;
                        cumm_delta_dn += (payoff(S_o_dn * exp(g * T + eps * sd)) +
                            payoff(S_o_dn * exp(g * T - eps * sd))) / 2;
                        cumm_rho_up += (payoff(S_o * exp(g_rho_up * 
                        T + eps * sd)) +
                            payoff(S_o * exp(g_rho_up * T - eps * sd))) / 2;
                        cumm_theta_up += (payoff(S_o * exp(g * (1 - incr) * T + eps * sigma * sqrt((1 - incr) * T))) +
                            payoff(S_o * exp(g * (1 - incr) * T - eps * sigma * sqrt((1 - incr) * T)))) / 2;
                        cumm_vega_up += (payoff(S_o * exp((r - q - sigma_up * sigma_up / 2) * T + eps * sigma_up * sqrt(T))) +
                            payoff(S_o * exp((r - q - sigma_up * sigma_up / 2) * T - eps * sigma_up * sqrt(T)))) / 2;
                    }
                }
                /* By linearity these greeks have been calculated once in terms of the accumulated
                   values. To estimate their stddevs, as we do for price, they would need to be moved inside
                   the loop over N_sims to and the squares of their values accumulated. The code as 
                   written is much more compact/readable...
                   */
                double delta_up = df * (cumm_delta_up - cumm_P) / (N_sims * incr * S_o);
                double delta_dn = df * (cumm_P - cumm_delta_dn) / (N_sims * incr * S_o);
                g_ptr->delta = (delta_dn + delta_up) / 2.0;
                g_ptr->gamma = (delta_up - delta_dn) / (incr * S_o);
                g_ptr->rho = df * (exp(-incr * r * T) * cumm_rho_up - cumm_P) / (N_sims * incr * r);
                g_ptr->theta = df * (exp(r * incr * T) * cumm_theta_up - cumm_P) / (N_sims * incr * T);
                g_ptr->vega = df * (cumm_vega_up - cumm_P) / (N_sims * incr * sigma);
            }
        } else { // Asian must simulate each time step:
            double dt = T / N_steps;
            double sd = sigma * sqrt(dt);
            double path_mean;
            double S;
            if (!g_ptr) { // nullptr -skip greeks
                if (!antithetic) {
                    for (int i = 0; i < N_sims; ++i) {
                        S = S_o;
                        path_mean = init_accum;
                        for (int j = 0; j < N_steps; ++j) {
                            double eps = norm(gen);
                            S *= exp(g * dt + eps * sd);
                            path_accum_fn(j + 1, S, path_mean);  // Could use (S_prev + S) / 2.
                        }
                        double P = (style == Terms::Style::AsianPrc) ? payoff(path_mean) : path_payoff(path_mean, S);
                        cumm_P += P;
                        ssq += P * P;
                    }
                } else { // Perform antithetic sampling:
                    // Introduce the antithetic path 'accumulators':
                    double path_mean_anti;
                    double S_anti;
                    for (int i = 1; i <= N_sims; ++i) {
                        // Initialize path accumulators prior to this(i) sim/run:
                        path_mean = init_accum;
                        path_mean_anti = init_accum;
                        S = S_o;
                        S_anti = S_o;
                        for (int j = 0; j < N_steps; ++j) {
                            double eps = norm(gen);
                            S *= exp(g * dt + eps * sd);
                            path_accum_fn(j + 1, S, path_mean);
                            S_anti *= exp(g * dt - eps * sd);
                            path_accum_fn(j + 1, S_anti, path_mean_anti);
                        }
                        double P = (style == Terms::Style::AsianPrc) ? (payoff(path_mean) + payoff(path_mean_anti)) / 2
                            : (path_payoff(path_mean, S) + path_payoff(path_mean_anti, S_anti)) / 2;
                        cumm_P += P;
                        ssq += P * P;
                    }
                }
            } else { // valid g_ptr, get greeks:
                double incr = 0.01;
                double sigma_up{ (1 + incr) * sigma };
                // Initialize accumulators before the N_sims:
                double cumm_delta_up{ 0 }, cumm_delta_dn{ 0 }, cumm_rho_up{ 0 },
                    cumm_theta_up{ 0 }, cumm_vega_up{ 0 };
                // S, path_mean and the following 7 path accumulators get reset for each of the N_sims by init_sim[&]():
                double S_rho_up, S_theta_up, S_vega_up;
                double path_mean_up, path_mean_dn;
                double pm_rho_up, pm_vega_up, pm_theta_up;
                if (!antithetic) {
                    std::function<void()> init_sim{ [&] () {
                        S = S_o; S_rho_up = S_o; S_theta_up = S_o; S_vega_up = S_o;
                        path_mean = init_accum; path_mean_up = init_accum; path_mean_dn = init_accum;
                        pm_rho_up = init_accum; pm_vega_up = init_accum;
                    } };
                    for (int i = 0; i < N_sims; ++i) {
                        init_sim();  // reset this sim run accumulators.
                        for (int j = 0; j < N_steps; ++j) {
                            double eps = norm(gen);
                            S *= exp(g * dt + eps * sd);
                            path_accum_fn(j + 1, S, path_mean);
                            path_accum_fn(j + 1, S * (1 + incr), path_mean_up);
                            path_accum_fn(j + 1, S * (1 - incr), path_mean_dn);
                            S_rho_up *= exp((g + incr * r) * dt + eps * sd);
                            path_accum_fn(j + 1, S_rho_up, pm_rho_up);
                            S_vega_up *= exp((r - q - sigma_up * sigma_up / 2) * dt + eps * sigma_up * sqrt(dt));
                            path_accum_fn(j + 1, S_vega_up, pm_vega_up);
                            if (j == N_steps - 2) { // special treatment -conditional has execution cost time(could put outside loop).
                                pm_theta_up = path_mean;
                                S_theta_up = S;
                            }
                        }
                        double P{};
                        if (style == Terms::Style::AsianPrc) {
                            P = payoff(path_mean);
                            cumm_delta_up += payoff(path_mean_up);
                            cumm_delta_dn += payoff(path_mean_dn);
                            cumm_rho_up += payoff(pm_rho_up);
                            cumm_vega_up += payoff(pm_vega_up);
                            cumm_theta_up += payoff(pm_theta_up);
                        } else {
                            P = path_payoff(path_mean, S);
                            cumm_delta_up += path_payoff(path_mean_up, S * (1 + incr));
                            cumm_delta_dn += path_payoff(path_mean_dn, S * (1 - incr));
                            cumm_rho_up += path_payoff(pm_rho_up, S);
                            cumm_vega_up += path_payoff(pm_vega_up, S);
                            cumm_theta_up += path_payoff(pm_theta_up, S_theta_up);
                        }
                        cumm_P += P;
                        ssq += P * P;
                    }
                } else {  // antithetic:
                    // Introduce the antithetic path 'accumulators':
                    double path_mean_anti;
                    double S_anti;
                    double S_rho_up_anti, S_theta_up_anti, S_vega_up_anti;
                    double path_mean_up_anti, path_mean_dn_anti;
                    double pm_rho_up_anti, pm_vega_up_anti, pm_theta_up_anti;

                    // init_sim_antith() called to initialize path accumulators prior to each of the N_sims:
                    std::function<void()> init_sim_antith{ [&] () {
                        S = S_o; S_rho_up = S_o; S_theta_up = S_o; S_vega_up = S_o;
                        path_mean = init_accum; path_mean_up = init_accum; path_mean_dn = init_accum;
                        pm_rho_up = init_accum; pm_vega_up = init_accum;
                        S_anti = S_o; S_rho_up_anti = S_o; S_theta_up_anti = S_o; S_vega_up_anti = S_o;
                        path_mean_anti = init_accum; path_mean_up_anti = init_accum; path_mean_dn_anti = init_accum;
                        pm_rho_up_anti = init_accum; pm_vega_up_anti = init_accum;
                    } };

                    for (int i = 0; i < N_sims; ++i) {
                        init_sim_antith();
                        for (int j = 0; j < N_steps; ++j) {
                            double eps = norm(gen);
                            S *= exp(g * dt + eps * sd);
                            path_accum_fn(j + 1, S, path_mean);
                            path_accum_fn(j + 1, S * (1 + incr), path_mean_up);
                            path_accum_fn(j + 1, S * (1 - incr), path_mean_dn);
                            S_rho_up *= exp((g + incr * r) * dt + eps * sd);
                            path_accum_fn(j + 1, S_rho_up, pm_rho_up);
                            S_vega_up *= exp((r - q - sigma_up * sigma_up / 2) * dt + eps * sigma_up * sqrt(dt));
                            path_accum_fn(j + 1, S_vega_up, pm_vega_up);
                            if (j == N_steps - 2) { // special treatment -conditional has execution cost time(could put outside loop).
                                pm_theta_up = path_mean;
                                S_theta_up = S;
                            }
                            // Same for antithetic path:
                            S_anti *= exp(g * dt - eps * sd);
                            path_accum_fn(j + 1, S_anti, path_mean_anti);
                            path_accum_fn(j + 1, S_anti * (1 + incr), path_mean_up_anti);
                            path_accum_fn(j + 1, S_anti * (1 - incr), path_mean_dn_anti);
                            S_rho_up_anti *= exp((g + incr * r) * dt - eps * sd);
                            path_accum_fn(j + 1, S_rho_up_anti, pm_rho_up_anti);
                            S_vega_up_anti *= exp((r - q - sigma_up * sigma_up / 2) * dt - eps * sigma_up * sqrt(dt));
                            path_accum_fn(j + 1, S_vega_up_anti, pm_vega_up_anti);
                            if (j == N_steps - 2) {
                                pm_theta_up_anti = path_mean_anti;
                                S_theta_up_anti = S_anti;
                            }
                        }
                        // Accumulate for this(i) sim/run:
                        double P{};
                        if (style == Terms::Style::AsianPrc) {
                            P = (payoff(path_mean) + payoff(path_mean_anti)) / 2;  // averaging antithetic path results.
                            cumm_delta_up += (payoff(path_mean_up) + payoff(path_mean_up_anti)) / 2;
                            cumm_delta_dn += (payoff(path_mean_dn) + payoff(path_mean_dn_anti)) / 2;
                            cumm_rho_up += (payoff(pm_rho_up) + payoff(pm_rho_up_anti)) / 2;
                            cumm_vega_up += (payoff(pm_vega_up) + payoff(pm_vega_up_anti)) / 2;
                            cumm_theta_up += (payoff(pm_theta_up) + payoff(pm_theta_up_anti)) / 2;
                        } else {
                            P = (path_payoff(path_mean, S) + path_payoff(path_mean_anti, S_anti)) / 2;
                            cumm_delta_up += (path_payoff(path_mean_up, S * (1 + incr)) + path_payoff(path_mean_up_anti, S_anti * (1 + incr))) / 2;
                            cumm_delta_dn += (path_payoff(path_mean_dn, S * (1 - incr)) + path_payoff(path_mean_dn_anti, S_anti * (1 - incr))) / 2;
                            cumm_rho_up += (path_payoff(pm_rho_up, S_rho_up) + path_payoff(pm_rho_up_anti, S_rho_up_anti)) / 2;
                            cumm_vega_up += (path_payoff(pm_vega_up, S_vega_up) + path_payoff(pm_vega_up_anti, S_vega_up_anti)) / 2;
                            cumm_theta_up += (path_payoff(pm_theta_up, S_theta_up) + path_payoff(pm_theta_up_anti, S_theta_up_anti)) / 2;
                        }
                        cumm_P += P;
                        ssq += P * P;
                    }

                }
                double delta_up = df * (cumm_delta_up - cumm_P) / (N_sims * incr * S_o);
                double delta_dn = df * (cumm_P - cumm_delta_dn) / (N_sims * incr * S_o);
                g_ptr->delta = (delta_dn + delta_up) / 2.0;
                g_ptr->gamma = (delta_up - delta_dn) / (incr * S_o);
                g_ptr->rho = df * (exp(-incr * r * T) * cumm_rho_up - cumm_P) / (N_sims * incr * r);
                g_ptr->theta = df * (exp(r * dt) * cumm_theta_up - cumm_P) / (N_sims * dt); // uses dt, not incr*T.
                g_ptr->vega = df * (cumm_vega_up - cumm_P) / (N_sims * incr * sigma);
            }
        }

        double mean_P = cumm_P / N_sims;
        sd_mean_PV = df * sqrt((ssq - N_sims * mean_P * mean_P)) / (N_sims - 1);
        double mean_PV = df * mean_P;

        /* // These had much higher variance than by using same sample paths for up/down estimates..
                if (g_ptr) {
                    double incr = 0.01;
                    double delta_up = (MC_prc_raw(N_sims, N_steps, r, style, T, payoff, path_payoff, path_accum_fn,
                        (1 + incr) * S_o, q, sigma, nullptr, antithetic, sd_mean_P) - mean_P) / (incr * S_o);
                    double delta_down = (mean_P - MC_prc_raw(N_sims, N_steps, r, style, T, payoff, path_payoff, path_accum_fn,
                        (1 - incr) * S_o, q, sigma, nullptr, antithetic, sd_mean_P)) / (incr * S_o);
                    g_ptr->delta = (delta_down + delta_up) / 2.0;
                    g_ptr->gamma = (delta_up - delta_down) / (incr * S_o);
                    g_ptr->rho = (MC_prc_raw(N_sims, N_steps, (1 + incr) * r, style, T, payoff, path_payoff, path_accum_fn,
                        S_o, q, sigma, nullptr, antithetic, sd_mean_P) - mean_P) / (incr * r);
                    g_ptr->vega = (MC_prc_raw(N_sims, N_steps, r, style, T, payoff, path_payoff, path_accum_fn,
                        S_o, q, (1 + incr) * sigma, nullptr, antithetic, sd_mean_P) - mean_P) / (incr * sigma);
                    g_ptr->theta = (MC_prc_raw(N_sims, N_steps, r, style, (1 - incr) * T, payoff, path_payoff, path_accum_fn,
                        S_o, q, sigma, nullptr, antithetic, sd_mean_P) - mean_P) / (incr * T);
                }
        */
        return mean_PV;
    }

    std::ostream& operator<< (std::ostream& os, const Terms& t) {
        switch (t.style()) {
        case Terms::Style::Euro: os << "European ";  break;
        case Terms::Style::Amer: os << "American ";  break;
        case Terms::Style::AsianPrc: os << "AsianPrc ";  break;
        case Terms::Style::AsianStrike: os << "AsianStrike ";  break;
        case Terms::Style::AsianExotic: os << "AsianExotic ";  break;
            // case Terms::Style::Other: os << "Other "; break;
        default: os.setstate(std::ios_base::failbit);
        }
        if (t.style() == Terms::Style::AsianPrc || t.style() == Terms::Style::AsianStrike
            || t.style() == Terms::Style::AsianExotic) {
            switch (t.avgType()) {
            case Terms::AvgType::Arith: os << "Arithmetic ";  break;
            case Terms::AvgType::Geom: os << "Geometric ";  break;
            case Terms::AvgType::Other: os << "Other ";  break;
            case Terms::AvgType::None: os << "None ";  break;
            default: os.setstate(std::ios_base::failbit);
            }
        }


        switch (t.type()) {
        case Terms::Type::Call: os << "Call  ";  break;
        case Terms::Type::Put: os << "Put   ";  break;
        case Terms::Type::Other: os << "Other   "; break;
        default: os.setstate(std::ios_base::failbit);
        }
        switch (t.style()) {
        case Terms::Style::Euro:
        case Terms::Style::Amer:
        case Terms::Style::AsianPrc: os << t.K();  break;
        default: os << '-';
        }
        os << ' ' << t.T();
        return os;
    }

    std::ostream& operator<< (std::ostream& os, const Asset& a) {
        os << a.S_o() << "  " << a.q() << "  " << a.sigma();
        return os;
    }

    std::ostream& operator<< (std::ostream& os, const Option& opt) {
        os << opt.terms() << ' ' << opt.asset();
        return os;
    }
}