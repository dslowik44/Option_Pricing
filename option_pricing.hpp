#ifndef OPTION_HDR
#define OPTION_HDR

#include <functional>
#include <iostream>
#include <random>

namespace OptionPricing {

    class Asset {
    public:
        Asset(double S_o, double q, double sigma);
        Asset(const Asset&) = default;

        double   S_o() const { return _S_o; }
        double     q() const { return _q; }
        double sigma() const { return _sigma; }

    private:
        double _S_o;
        double _q;
        double _sigma;
    };


    class Terms {
    public:
        enum class Style { Euro, Amer, AsianPrc, AsianStrike, AsianExotic };
        enum class Type { Call, Put, Other };
        enum class AvgType { Arith, Geom, Other, None };

        static std::function<void(int, double, double&)> arith_path_accum_fn;
        static std::function<void(int, double, double&)> geom_path_accum_fn;

        Terms(Style, Type, double T, double K, AvgType = AvgType::Arith);
        Terms(Style, Type, double T, AvgType = AvgType::Arith);
        Terms(Style, Type, double T, std::function<double(double)> payoff);
        Terms(Style, Type, double T, std::function<double(double, double)> path_payoff,
            std::function<void(int, double, double&)> path_accum, double init_accum);
        Terms(const Terms&) = default;
        Terms() = delete;

        Style style() const { return _style; }
        Type type() const { return _type; }
        double K() const { return _K; }
        double T() const { return _T; }
        AvgType avgType() const { return _avg_type; }
        std::function<double(double)> payoff_fn() const { return _payoff; }
        std::function<double(double, double)> path_payoff() const { return _path_payoff; }
        std::function<void(int, double, double&)> path_accum() const { return _path_accum_fn; }
        double init_accum() const { return _init_accum; }

    private:
        const Style _style;
        const Type _type;
        const double _T;
        const double _K;
        AvgType _avg_type;
        std::function<double(double)> _payoff{};
        std::function<double(double, double)> _path_payoff{};
        std::function<void(int, double, double&)> _path_accum_fn{};
        double _init_accum{ 0.0 };
    };


    class Option {
    public:
        Option(const Terms&, const Asset&);

        double btree_prc(int N_steps, double r, bool get_greeks);
        double MC_prc(int N_sims, int N_steps, double r, bool get_greeks,
            bool antithetic, double& sd_prc);

        struct Greeks {
            double delta{ 0 };
            double gamma{ 0 };
            double theta{ 0 };
            double   rho{ 0 };
            double  vega{ 0 };
        };

        Terms  terms() const { return _terms; }
        Asset  asset() const { return _asset; }
        double price() const { return _prc; }
        Greeks all_greeks() const { return _greeks; }
        double delta() const { return _greeks.delta; }
        double gamma() const { return _greeks.gamma; }
        double theta() const { return _greeks.theta; }
        double   rho() const { return _greeks.rho; }
        double  vega() const { return _greeks.vega; }

    private:
        const Terms _terms;
        const Asset _asset;
        double _prc;
        Greeks _greeks;
    };

    double btree_prc_opt(Option&, int N_steps, double r, bool greeks);

    double btree_prc_terms_asset(const Terms&, const Asset&, int N_steps, double r,
        Option::Greeks*); // Note user danger here.

    // All binary tree pricing calls lead here:
    double btree_prc_raw(int N_steps, double r,
        Terms::Style style, double T, std::function<double(double)>,
        double S_o, double q, double sigma, Option::Greeks*);


    double MC_prc_opt(Option&, int N_sims, int N_steps, double r, bool greeks,
        bool antithetic, double& sd_prc);

    double MC_prc_terms_asset(const Terms&, const Asset&, int N_sims, int N_steps, double r,
        Option::Greeks*, bool antithetic, double& sd_prc);

    // All MC pricing calls lead here:
    double MC_prc_raw(int N_sims, int N_steps, double r,
        Terms::Style style, double T, const std::function<double(double)>& payoff,
        const std::function<double(double, double)>& path_payoff,
        const std::function<void(int, double, double&)>& path_accum, double init_accum,
        double S_o, double q, double sigma, Option::Greeks*,
        bool antithetic, double& sd_prc);

    namespace MC_InternalUtilities {              
        double MC_prc_Euro(int, double, double, const std::function<double(double)>&,
            double, double, double, Option::Greeks*, bool, double&);

        double MC_prc_Asian(int, const int, double, Terms::Style, double, const std::function<double(double)>&,
                const std::function<double(double, double)>&, const std::function<void(int, double, double&)>&, double,
                double, double, double, Option::Greeks*, bool, double&);

        std::mt19937& MC_get_rndm_gen();
    }

    std::ostream& operator<< (std::ostream&, const Asset&);
    std::ostream& operator<< (std::ostream&, const Terms&);
    std::ostream& operator<< (std::ostream&, const Option&);

}

#endif
