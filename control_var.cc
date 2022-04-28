#include "option_pricing.hpp"
#include <iostream>
#include <cmath>
#include <valarray>
#include <iomanip>
#include <fstream>

using namespace OptionPricing;

//using Terms::Style::Amer, Terms::Type::Put;
using std::cout;
using std::endl;
using std::valarray;


int main() {
    // Construct the Amer and Euro Put options for a particular parameter set:
    double S_o = 40.0, r = 0.0488, q = 0.0;
    double K = 40.0, T = 7 / 12., sigma = 0.2;
    Asset underlying{ S_o, q, sigma };
    Terms terms_amer_put{ Terms::Style::Amer, Terms::Type::Put, T, K };
    Terms terms_euro_put{ Terms::Style::Euro, Terms::Type::Put, T, K };
    Option amer_put(terms_amer_put, underlying);
    Option euro_put(terms_euro_put, underlying);

    // Get their prices on the training set of binomial trees having N_steps on [20, 100] into
    // the arrays amer_prcs, euro_prcs:
    bool get_greeks = false;
    constexpr int start_N = 20, end_N = 100;
    constexpr int N_points{ end_N - start_N + 1 }; // Number of training points.

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
    cout << "Training Set (binomial trees having N_steps on [" << start_N << " : "
        << end_N << "]) Results:" << endl;
    cout << "Optimal c: " << c << endl;
    cout << "Corr(Amer, Euro): " << cov_a_e / (sqrt(var_amer) * sqrt(var_euro)) << endl;

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
    cout << "         c                  Var | stddev" << endl;
    cout << "None                 " << var_amer << '|' << sqrt(var_amer) << endl;
    cout << "Hull-White, c=1      " << var_amer_HW << '|' << sqrt(var_amer_HW) << endl;
    cout << "Optimal c=Cov/Var    " << var_amer_opt << '|' << sqrt(var_amer_opt) << endl;
    cout << "Mean of raw American prices: " << amer_mean << endl;
    cout << "Mean of Hull-White prices:   " << amer_mean_HW << endl;
    cout << "Mean of opt c* prices:       " << amer_mean_opt << endl;
    cout << "Mean of raw European prices: " << euro_mean << endl;
    cout << "Black-Scholes European price: " << BS_prc << endl;


    // Out of sample test on binomial trees same number, N_points, of points as training sample
    // but all trees having N_steps >> end_N.
    // Test grid: [1000 : 1000 + N_points]
    int test_start_N = 4000;
    int test_end_N = 4000 + N_points - 1;
    cout << "\nTest set (binomial trees having N_steps on [" << test_start_N << " : "
        << test_end_N << "]) Results:" << endl;

    for (int N_steps = test_start_N; N_steps < test_end_N + 1; ++N_steps) {
        amer_prcs[N_steps - test_start_N] = amer_put.btree_prc(N_steps, r, get_greeks);
        euro_prcs[N_steps - test_start_N] = euro_put.btree_prc(N_steps, r, get_greeks);
    }

    // Calculate their sample means, vars and covariance:
    amer_mean = amer_prcs.sum() / N_points;
    euro_mean = euro_prcs.sum() / N_points;

    var_amer = ((amer_prcs - amer_mean) * (amer_prcs - amer_mean)).sum() / (N_points - 1);
    var_euro = ((euro_prcs - euro_mean) * (euro_prcs - euro_mean)).sum() / (N_points - 1);
    cov_a_e = ((amer_prcs - amer_mean) * (euro_prcs - euro_mean)).sum() / (N_points - 1);
    // Get the theoretical optimum mix of the Euro option control variate (for this parameter set):
    double c_testset = cov_a_e / var_euro;
    cout << "Optimal c (in test set): " << c_testset << endl;
    cout << "Test set Corr(Amer, Euro): " << cov_a_e / (sqrt(var_amer) * sqrt(var_euro)) << endl;

    // Get the variance of Amer option prices via the optimal, and the Hull, White (c==1) control
    // variate admixture:
    amer_prcs_HW = amer_prcs - (euro_prcs - BS_prc);
    amer_prcs_opt = amer_prcs - c * (euro_prcs - BS_prc); // Use c from training set!
    amer_mean_HW = amer_prcs_HW.sum() / N_points;
    amer_mean_opt = amer_prcs_opt.sum() / N_points;
    var_amer_HW = ((amer_prcs_HW - amer_mean_HW) * (amer_prcs_HW - amer_mean_HW)).sum() / (N_points - 1);
    var_amer_opt = ((amer_prcs_opt - amer_mean_opt) * (amer_prcs_opt - amer_mean_opt)).sum() / (N_points - 1);

    std::fstream outfile("prcs_4000-4080.csv", std::ios_base::out);
    for (int N_steps = test_start_N; N_steps < test_end_N + 1; ++N_steps) {
        outfile << N_steps << "," << amer_prcs[N_steps - test_start_N] << ',' 
                << euro_prcs[N_steps - test_start_N] << ',' 
                << std::setprecision(10) << amer_prcs_opt[N_steps - test_start_N] << endl;
    }

    cout << "Var and stddev of Z = A -c(E - BS_prc):" << endl;
    cout << "         c                  Var | stddev" << endl;
    cout << "None                 " << var_amer << '|' << sqrt(var_amer) << endl;
    cout << "Hull-White, c=1      " << var_amer_HW << '|' << sqrt(var_amer_HW) << endl;
    cout << "Optimal c=Cov/Var    " << var_amer_opt << '|' << sqrt(var_amer_opt) << endl;
    cout << "Mean of raw American prices: " << amer_mean << endl;
    cout << "Mean of Hull-White prices:   " << amer_mean_HW << endl;
    cout << "Mean of opt c* prices:       " << amer_mean_opt << endl;
    cout << "Mean of raw European prices: " << euro_mean << endl;
   
   outfile.close();

   // write Amer, Euro and Optimal prices for N_steps on [1000, 20000]:
   outfile.open("prcs_1000-20000.csv", std::ios_base::out);
   for (int N_steps = 1000; N_steps < 20002; ++N_steps) {
       if (N_steps % 1000 < 2) {
        double a = amer_put.btree_prc(N_steps, r, get_greeks);
        double e = euro_put.btree_prc(N_steps, r, get_greeks);
        outfile << N_steps << "," << a << ',' << e << ',' << a - c * (e - BS_prc) << endl;
        // cout << N_steps << ' ' << std::setprecision(10) << a - c * (e - BS_prc) << endl;
        // cout << N_steps << ' ' << std::setprecision(10) << a << endl;
       }
    }
   

}