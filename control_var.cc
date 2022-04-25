#include "option_pricing.hpp"
// #include "fexcpt.hpp"
#include <iostream>
#include <cmath>
#include <memory>
#include <numeric>
#include <iterator>
#include <valarray>

using namespace OptionPricing;

//using Terms::Style::Amer, Terms::Type::Put; 
using std::cout;
using std::endl;
using std::valarray;


int main() {
    // Construct the Amer and Euro Put options:
    double S_o = 40.0, r = 0.0488, q = 0.0;
    double K = 40.0, T = 7 / 12., sigma = 0.2;
    Asset underlying{ S_o, q, sigma };
    Terms terms_amer_put{ Terms::Style::Amer, Terms::Type::Put, T, K };
    Terms terms_euro_put{ Terms::Style::Euro, Terms::Type::Put, T, K };
    Option amer_put(terms_amer_put, underlying);
    Option euro_put(terms_euro_put, underlying);
    
    // Get their prices calculated on binomial trees having N_steps on [20, 100] into 
    // the arrays amer_prcs, euro_prcs:
    bool get_greeks = false;
    constexpr int start_N = 20, end_N = 100;
    constexpr int N_points { end_N - start_N + 1 };

    valarray<double> amer_prcs(N_points);
    valarray<double> euro_prcs(N_points);
    
    
    for (int N_steps = start_N; N_steps < end_N + 1; ++N_steps) {
        amer_prcs[N_steps - start_N] = amer_put.btree_prc(N_steps, r, get_greeks);
        euro_prcs[N_steps - start_N] = euro_put.btree_prc(N_steps, r, get_greeks);
    }

    // Calculate their sample means, vars and covariance:
    double amer_mean = amer_prcs.sum() / N_points;
    double euro_mean = euro_prcs.sum() / N_points;

    double var_amer = ((amer_prcs - amer_mean) * (amer_prcs - amer_mean)).sum() / (N_points - 1);
    double var_euro = ((euro_prcs - euro_mean) * (euro_prcs - euro_mean)).sum() / (N_points - 1);
    double cov_a_e = ((amer_prcs - amer_mean) * (euro_prcs - euro_mean)).sum() / (N_points - 1);


    // Get the theoretical optimum mix of the Euro option control variate (for this parameter set):
    double c = cov_a_e / var_euro;

    // Get the Black-Scholes Euro option price:
    double d_1 = (std::log(S_o / K) + (r - q + sigma * sigma / 2) * T) / (sigma * std::sqrt(T));
    double d_2 = d_1 - sigma * std::sqrt(T);
    double BS_prc = std::exp(-r * T) * K * std::erfc(d_2 / std::sqrt(2)) / 2
                    - S_o * std::erfc(d_1 / std::sqrt(2)) / 2;

    // Get the variance of Amer option prices via the optimal, and the Hull, White (c==1) control
    // variate admixture:
    valarray<double> amer_prcs_HW = amer_prcs - (euro_prcs - BS_prc);
    valarray<double> amer_prcs_opt = amer_prcs - c * (euro_prcs - BS_prc);
    double amer_mean_HW = amer_prcs_HW.sum() / N_points;
    double amer_mean_opt = amer_prcs_opt.sum() / N_points;
    double var_amer_HW = ((amer_prcs_HW - amer_mean_HW) * (amer_prcs_HW - amer_mean_HW)).sum() / (N_points - 1);
    double var_amer_opt = ((amer_prcs_opt - amer_mean_opt) * (amer_prcs_opt - amer_mean_opt)).sum() / (N_points - 1);
    cout << "Control Variate       Var | stddev" << endl;
    cout << "None             " << var_amer << '|' << sqrt(var_amer) << endl;
    cout << "Hull-White, c=1  " << var_amer_HW << '|' << sqrt(var_amer_HW) << endl;
    cout << "Optimal c=Cov/Var" << var_amer_opt << '|' << sqrt(var_amer_opt) << endl;



}

