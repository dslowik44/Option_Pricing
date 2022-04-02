#include "option_pricing.hpp"
// #include "fexcpt.hpp"
#include <iostream>
#include <cmath>


using namespace OptionPricing;

//using Terms::Style::Amer, Terms::Type::Put; 
using std::cout;
using std::endl;

double straddle_payoff(double S) {
    if (S > 40.2)
        return S - 40.2;
    else
        return 40.2 - S;
}

int main() {
    // clear_all_fe();

    // Example 21.1 Hull. First create the American put option from its terms and underlying asset:
    double K = 50.0, T = 5 / 12.;
    Terms terms{ Terms::Style::Amer, Terms::Type::Put, T, K };
    double S_o = 50.0, q = 0.0, sigma = 0.40;
    Asset underlying{ S_o, q, sigma };
    Option amer_put(terms, underlying);

    // Now price the option on a binary tree:
    int N_steps = 5;
    double r = 0.10;
    bool get_greeks = false;
    amer_put.btree_prc(N_steps, r, get_greeks);
    // flag_fe("After calling amer_put.btree_prc", false);
    cout << "Hull, Example 21.1:\n"
        << "Binary tree pricing of option:\n"
        << "r   Style    Type  K   T     S_o q sigma\n"
        << r << ' ' << amer_put << '\n'
        << "N  price:\n"
        << N_steps << ' ' << amer_put.price() << endl;
    for (int N_steps : {30, 50, 100, 500}) {
        amer_put.btree_prc(N_steps, r, get_greeks);
        cout << N_steps << ' ' << amer_put.price() << endl;
    }

    // Replicate Fig 21.4: price fluctuations and convergence as N_steps increases:
    constexpr int NN{ 51 };
    double prcs[NN]{ 0 };
    cout << "\n Fig. 21.4, N vs Price:\nN  Price:\n";
    for (int N_steps = 2; N_steps < NN; ++N_steps) {
        prcs[N_steps] = amer_put.btree_prc(N_steps, r, get_greeks);
        cout << N_steps << ' ' << prcs[N_steps] << endl;
    }

    // Replicate Example 21.2, Estimates of delta, gamma, theta, rho and vega on amer_put:
    N_steps = 50;
    get_greeks = true;
    amer_put.btree_prc(N_steps, r, get_greeks);
    //flag_fe("After calling amer_put.btree_prc final time", false);
    cout << "\n Hull, Example 21.2, Greeks found with binary pricing tree with N_steps: " << N_steps << '\n'
        << "r   Style    Type  K   T     S_o q sigma\n"
        << r << ' ' << amer_put << ":\n"
        << "price: " << amer_put.price() << '\n'
        << "delta: " << amer_put.delta() << '\n' << "gamma: " << amer_put.gamma() << '\n'
        << "theta: " << amer_put.theta() << '\n' << "rho: " << amer_put.rho() << '\n'
        << "vega: " << amer_put.vega() << endl;

    // Cash-or-nothing call, see https://en.wikipedia.org/wiki/Binary_option
    std::function<double(double)> payoff_cash_call;
    payoff_cash_call = [K](double s) {
        if (s >= K)
            return 1.0;
        else
            return 0.0;
    };
    Terms terms_cash_call{ Terms::Style::Euro, Terms::Type::Other, T, payoff_cash_call };
    Option euro_cash_call(terms_cash_call, underlying);
    N_steps = 1000;
    euro_cash_call.btree_prc(N_steps, r, get_greeks);
    cout << "\n European cash-or-nothing call with N_steps: " << N_steps << '\n'
        << "r   Style    Type     K   T     S_o q sigma\n"
        << r << ' ' << euro_cash_call << ":\n"
        << "price: " << euro_cash_call.price() << '\n'
        << "delta: " << euro_cash_call.delta() << '\n' << "gamma: " << euro_cash_call.gamma() << '\n'
        << "theta: " << euro_cash_call.theta() << '\n' << "rho: " << euro_cash_call.rho() << '\n'
        << "vega: " << euro_cash_call.vega() << '\n';

    double d_1 = (std::log(S_o / K) + (r - q + sigma * sigma / 2) * T) / (sigma * std::sqrt(T));
    double d_2 = d_1 - sigma * std::sqrt(T);
    double BS_price = std::exp(-r * T) * std::erfc(-d_2 / std::sqrt(2)) / 2;
    cout << "Black-Scholes-Merton price: " << BS_price << endl;

    // Asset-or-nothing call:
    std::function<double(double)> payoff_asset_call;
    payoff_asset_call = [K](double s) {
        if (s > K)
            return s;
        else
            return 0.0;
    };
    Terms terms_asset_call{ Terms::Style::Euro, Terms::Type::Other, T, payoff_asset_call };
    Option euro_asset_call(terms_asset_call, underlying);
    N_steps = 10000;
    euro_asset_call.btree_prc(N_steps, r, get_greeks);
    cout << "\n European asset-or-nothing call with N_steps: " << N_steps << '\n'
        << "r   Style    Type     K   T     S_o q sigma\n"
        << r << ' ' << euro_asset_call << ":\n"
        << "price: " << euro_asset_call.price() << '\n'
        << "delta: " << euro_asset_call.delta() << '\n' << "gamma: " << euro_asset_call.gamma() << '\n'
        << "theta: " << euro_asset_call.theta() << '\n' << "rho: " << euro_asset_call.rho() << '\n'
        << "vega: " << euro_asset_call.vega() << '\n';

    BS_price = S_o * std::exp(-q * T) * std::erfc(-d_1 / std::sqrt(2)) / 2;
    cout << "Black-Scholes-Merton price: " << BS_price << endl;


    // Euro option via binary tree:
    T = 3 / 12.0;
    K = 43.0;
    Terms terms_euro{ Terms::Style::Euro, Terms::Type::Put, T, K };
    S_o = 40.2;
    q = 0.2;
    sigma = 0.32;
    Asset underlying_euro{ S_o, q, sigma };
    Option euro_put(terms_euro, underlying_euro);
    r = 0.07;
    N_steps = 1000;
    get_greeks = true;
    euro_put.btree_prc(N_steps, r, get_greeks);
    cout << "\nOn Binary pricing tree with N_steps: " << N_steps << '\n'
        << "Style    Class  K   T     S_o q sigma\n"
        << euro_put << ":\n"
        << "price: " << euro_put.price() << '\n'
        << "delta: " << euro_put.delta() << '\n' << "gamma: " << euro_put.gamma() << '\n'
        << "theta: " << euro_put.theta() << '\n' << "rho: " << euro_put.rho() << '\n'
        << "vega: " << euro_put.vega() << endl;

    // Same Euro option via Monte Carlo pricer:
    double sd_prc{ 0.0 };
    int N_sims{ 10000 };
    bool antithetic{ true };
    get_greeks = true;
    euro_put.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
    cout << "\nEuro put with Monte Carlo pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
        << "Style    Class  K   T     S_o q sigma\n"
        << euro_put << ":\n"
        << "price: " << euro_put.price() << ", std_dev: " << sd_prc << '\n'
        << "delta: " << euro_put.delta() << '\n' << "gamma: " << euro_put.gamma() << '\n'
        << "theta: " << euro_put.theta() << '\n' << "rho: " << euro_put.rho() << '\n'
        << "vega: " << euro_put.vega() << endl;

    // At the money Straddle option via Monte Carlo pricer:
    Terms terms_straddle{ Terms::Style::Euro, Terms::Type::Other, T, straddle_payoff };
    Option euro_straddle(terms_straddle, underlying_euro);
    N_sims = 10000;
    antithetic = false;  // setting to true does not decrease std_dev of price found.
    get_greeks = true;
    sd_prc = 0.0;
    euro_straddle.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
    cout << "\nEuro Straddle with MC pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
        << "r    Style    Class    K T    S_o  q  sigma\n"
        << r << ' ' << euro_straddle << ":\n"
        << "price: " << euro_straddle.price() << ", std_dev: " << sd_prc << '\n'
        << "delta: " << euro_straddle.delta() << '\n' << "gamma: " << euro_straddle.gamma() << '\n'
        << "vega: " << euro_straddle.vega() << '\n' << "theta: " << euro_straddle.theta() << '\n'
        << "rho: " << euro_straddle.rho() << endl;

    // Deep in the money call where antithetic sampling decreases price std_dev by 9x:
    T = 3 / 12.0;
    K = 40.0;
    Terms terms_euro_call{ Terms::Style::Euro, Terms::Type::Call, T, K };
    S_o = 60;
    q = 0.2;
    sigma = 0.32;
    Asset underlying_euro_call{ S_o, q, sigma };
    Option euro_call(terms_euro_call, underlying_euro_call);
    r = 0.07;
    N_steps = 1000;
    get_greeks = true;
    N_sims = 100000;
    antithetic = true;
    sd_prc = 0.0;
    euro_call.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
    cout << "\nEuro call with MC pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
        << "r    Style    Class  K T   S_o  q  sigma\n"
        << r << ' ' << euro_call << ":\n"
        << "price: " << euro_call.price() << ", std_dev: " << sd_prc << '\n'
        << "delta: " << euro_call.delta() << '\n' << "gamma: " << euro_call.gamma() << '\n'
        << "vega: " << euro_call.vega() << '\n' << "theta: " << euro_call.theta() << '\n'
        << "rho: " << euro_call.rho() << endl;

    // AsianPrc Put default AvgType::Arith
    T = 6 / 12.0;
    K = 50.0;
    Terms terms_asianprc_put{ Terms::Style::AsianPrc, Terms::Type::Put, T, K };
    S_o = 49.0;
    q = 0.0;
    sigma = 0.3;
    Asset underlying_asian{ S_o, q, sigma };
    Option asianprc_put(terms_asianprc_put, underlying_asian);
    r = 0.07;
    N_steps = 1000;
    N_sims = 10000;
    antithetic = true;
    get_greeks = true;
    sd_prc = 0.0;
    asianprc_put.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
    cout << "\nAsianPrc put with Monte Carlo pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
        << "r    Style    AvgType   Class  K  T  S_o  q   sigma\n"
        << r << ' ' << asianprc_put << ":\n"
        << "price: " << asianprc_put.price() << ", std_dev: " << sd_prc << '\n'
        << "delta: " << asianprc_put.delta() << '\n' << "gamma: " << asianprc_put.gamma() << '\n'
        << "theta: " << asianprc_put.theta() << '\n' << "rho: " << asianprc_put.rho() << '\n'
        << "vega: " << asianprc_put.vega() << endl;

    // Price an AsianPrc Put option with Geometric path average price:
    Terms terms_asianprc_put_geom{ Terms::Style::AsianPrc, Terms::Type::Put, T, K, Terms::AvgType::Geom };
    Option asianprc_put_geom(terms_asianprc_put_geom, underlying_asian);
    N_steps = 1000;
    N_sims = 10000;
    antithetic = true;
    get_greeks = true;
    sd_prc = 0.0;
    asianprc_put_geom.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
    cout << "\nAsianPrc put with Monte Carlo pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
        << "r    Style        AvgType   Class  K  T  S_o  q  sigma\n"
        << r << ' ' << asianprc_put_geom << ":\n"
        << "price: " << asianprc_put_geom.price() << ", std_dev: " << sd_prc << '\n'
        << "delta: " << asianprc_put_geom.delta() << '\n' << "gamma: " << asianprc_put_geom.gamma() << '\n'
        << "theta: " << asianprc_put_geom.theta() << '\n' << "rho: " << asianprc_put_geom.rho() << '\n'
        << "vega: " << asianprc_put_geom.vega() << endl;

    // next:

        std::function<double(double, double)> path_payoff{
            [](double path_avg, double) { return path_avg; }
        };

        S_o = 49.0;
        std::function<void(int, double, double&)> path_accum_tanh{
            [S_o](int idx, double S, double& accum) { accum += (std::tanh(S / S_o) - accum) / idx; }
        };

        T = 6 / 12.0;
        Terms terms_asian_tanh_call{ Terms::Style::AsianExotic, Terms::Type::Other, T, path_payoff,
                                     path_accum_tanh, 0.0 };
        q = 0.0;
        sigma = 0.3;
        Asset underlying_asian2{ S_o, q, sigma };
        Option asian_tanh_call(terms_asian_tanh_call, underlying_asian2);
        r = 0.07;
        N_steps = 1000;
        N_sims = 10000;
        antithetic = true;
        get_greeks = true;
        sd_prc = 0.0;
        asian_tanh_call.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
        cout << "\nAsianExotic tanh_path_avg 0_call MC pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
            << "r    Style     AvgType Class  K  T  S_o  q  sigma\n"
            << r << ' ' << asian_tanh_call << ":\n"
            << "price: " << asian_tanh_call.price() << ", std_dev: " << sd_prc << '\n'
            << "delta: " << asian_tanh_call.delta() << '\n' << "gamma: " << asian_tanh_call.gamma() << '\n'
            << "theta: " << asian_tanh_call.theta() << '\n' << "rho: " << asian_tanh_call.rho() << '\n'
            << "vega: " << asian_tanh_call.vega() << endl;


        std::function<double(double, double)> path_payoff_sin_tanh_abs{
            [](double path_avg, double S) { return std::abs(std::sin((path_avg - std::tanh(S)) / (std::tanh(S) / 10.0))); }
        };

        Terms terms_asian_tanh_sin_abs{ Terms::Style::AsianExotic, Terms::Type::Other, T, path_payoff_sin_tanh_abs,
                                     path_accum_tanh, 0.0 };
        Option asian_tanh_sin_abs{ terms_asian_tanh_sin_abs, underlying_asian2 };
        antithetic = true;
        get_greeks = true;
        sd_prc = 0.0;
        asian_tanh_sin_abs.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
        cout << "\nAsianExotic abs(sin(tanh_path_avg)) MC pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
            << "r    Style     AvgType Class  K  T  S_o  q  sigma\n"
            << r << ' ' << asian_tanh_sin_abs << ":\n"
            << "price: " << asian_tanh_sin_abs.price() << ", std_dev: " << sd_prc << '\n'
            << "delta: " << asian_tanh_sin_abs.delta() << '\n' << "gamma: " << asian_tanh_sin_abs.gamma() << '\n'
            << "theta: " << asian_tanh_sin_abs.theta() << '\n' << "rho: " << asian_tanh_sin_abs.rho() << '\n'
            << "vega: " << asian_tanh_sin_abs.vega() << endl;


        std::function<double(double, double)> path_payoff_sin_tanh{
            [](double path_avg, double S) { return std::sin((path_avg - std::tanh(S)) / (std::tanh(S) / 10.0)); }
        };

        Terms terms_asian_tanh_sin{ Terms::Style::AsianExotic, Terms::Type::Other, T, path_payoff_sin_tanh,
                                     path_accum_tanh, 0.0 };
        Option asian_tanh_sin{ terms_asian_tanh_sin, underlying_asian2 };
        antithetic = true;
        get_greeks = true;
        sd_prc = 0.0;
        asian_tanh_sin.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
        cout << "\nAsianExotic sin(tanh_path_avg) MC pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
            << "r    Style     AvgType Class  K  T  S_o  q  sigma\n"
            << r << ' ' << asian_tanh_sin << ":\n"
            << "price: " << asian_tanh_sin.price() << ", std_dev: " << sd_prc << '\n'
            << "delta: " << asian_tanh_sin.delta() << '\n' << "gamma: " << asian_tanh_sin.gamma() << '\n'
            << "theta: " << asian_tanh_sin.theta() << '\n' << "rho: " << asian_tanh_sin.rho() << '\n'
            << "vega: " << asian_tanh_sin.vega() << endl;

        /* need euro_put -> euro_straddle:
            euro_put.btree_prc(N_steps, r, get_greeks);
            cout << "\nOn Binary pricing tree with N_steps: " << N_steps << '\n'
                << "Style    Class  K   T     S_o q sigma\n"
                << euro_put << ":\n"
                << "price: " << euro_put.price() << '\n'
                << "delta: " << euro_put.delta() << '\n' << "gamma: " << euro_put.gamma() << '\n'
                << "theta: " << euro_put.theta() << '\n' << "rho: " << euro_put.rho() << '\n'
                << "vega: " << euro_put.vega() << endl;

            sd_prc = 0.0;
            euro_put.MC_prc(N_sims, N_steps, r, get_greeks, antithetic, sd_prc);
            cout << "\nEuro put with Monte Carlo pricer with N_sims: " << N_sims << ", N_steps: " << N_steps << '\n'
                << "Style    Class  K   T     S_o q sigma\n"
                << euro_put << ":\n"
                << "price: " << euro_put.price() << ", std_dev: " << sd_prc << '\n'
                << "delta: " << euro_put.delta() << '\n' << "gamma: " << euro_put.gamma() << '\n'
                << "theta: " << euro_put.theta() << '\n' << "rho: " << euro_put.rho() << '\n'
                << "vega: " << euro_put.vega() << endl;
        */
}

/*
    // Ex 21.3:
    K = 300.0;
    T = 4. / 12.;
    Terms terms2{ Terms::Style::Amer, Terms::Type::Call, T, K };
    S_o = 300.0;
    r = 0.08;
    q = r;
    sigma = 0.30;
    Asset idx_fut{ S_o, q, sigma };
    Option amer_call_idx_fut{ terms2, idx_fut };
    N_steps = 1000;
    amer_call_idx_fut.btree_prc(N_steps, r, true);
    cout << "\n Hull, Example 21.3, American Call on index futures contract, N_steps: " << N_steps << '\n'
        << "r   Style    Type  K   T     S_o q sigma\n"
        << r << ' ' << amer_call_idx_fut << ":\n"
        << "price: " << amer_call_idx_fut.price() << '\n'
        << "delta: " << amer_call_idx_fut.delta() << '\n' << "gamma: " << amer_call_idx_fut.gamma() << '\n'
        << "theta: " << amer_call_idx_fut.theta() << '\n' << "rho: " << amer_call_idx_fut.rho() << '\n'
        << "vega: " << amer_call_idx_fut.vega() << endl;

    // Call via namespace functions:
    cout << "\n Called via namespace function btree_prc_opt: \n"
        << "price: " << btree_prc_opt(amer_call_idx_fut, N_steps, r, true) << '\n'
        << "delta: " << amer_call_idx_fut.delta() << '\n' << "gamma: " << amer_call_idx_fut.gamma() << '\n'
        << "theta: " << amer_call_idx_fut.theta() << '\n' << "rho: " << amer_call_idx_fut.rho() << '\n'
        << "vega: " << amer_call_idx_fut.vega() << endl;

    Option::Greeks greeks{};
    cout << "\n Called via namespace function btree_prc_terms_assets: \n"
        << "price: " << btree_prc_terms_assets(terms2, idx_fut, N_steps, r, &greeks) << '\n'
        << "delta: " << greeks.delta << '\n' << "gamma: " << greeks.gamma << '\n'
        << "theta: " << greeks.theta << '\n' << "rho: " << greeks.rho << '\n'
        << "vega: " << greeks.vega << endl;

    cout << "\n Called via namespace function btree_prc_raw: \n"
        << "price: " << btree_prc_raw(N_steps, r, terms2.style(), T, terms2.payoff_fn(), S_o, q, sigma, &greeks) << '\n'
        << "delta: " << greeks.delta << '\n' << "gamma: " << greeks.gamma << '\n'
        << "theta: " << greeks.theta << '\n' << "rho: " << greeks.rho << '\n'
        << "vega: " << greeks.vega << endl;

*/
