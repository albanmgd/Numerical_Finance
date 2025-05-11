#include <iostream>
#include <vector>

#include "../Pricer/BermudeanBasketOption.h"
#include "KakutaniSequence.h"
#include "../Pricer/EuropeanBasketOption.h"
#include "../SDE/BlackScholesND.h"
#include "../RandomGenerator/LinearCongruential.h"
#include "../RandomGenerator/EcuyerCombined.h"
#include "../RandomGenerator/FiniteSet.h"
#include "../RandomGenerator/Exponential.h"
#include "../RandomGenerator/Normal.h"
#include "../RandomGenerator/Poisson.h"
#include "../RandomGenerator/PAdic.h"
#include "../PDE/PDEGrid2D.h"
#include "../SDE/BlackScholes1D.h"
#include "../SDE/BSEuler1D.h"
#include "../SDE/BrownianND.h"
#include "../SDE/BSEulerND.h"
#include "../SDE/BlackScholesND.h"
#include "../Utils/CSVWriter.h"

void EuroBasket();
void BermudBasket();
void BasketKakutani();

int main() {
//    EuroBasket();
//      BermudBasket();
    BasketKakutani();
}

void BermudBasket() {
    // first we set the parameters
    int nb_assets = 3;
    double maturity = 1.0;
    double strike = 60;
    size_t nb_steps = 365;
    // size_t nb_sim = 1e4;
    vector<double> spots = {100,50,60};
    vector<double> vols = {0.10,0.25,0.16};
    double rate = 0.05;
    vector<double> weights = {0.10, 0.7, 0.2};
    vector<vector<double>> correl_mat(nb_assets, vector<double>(nb_assets, 0.1));
    for (int i = 0; i < nb_assets; i++) {
        correl_mat[i][i] = 1.0;
    }
    bool use_control_variate = true;
    bool use_antithetic = false;

    // we start the pricing
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0.,1., Unif);
    BermudeanBasketOption bermud_basket_opt(nb_assets, strike, maturity, rate, spots, vols, weights, correl_mat, NormBox, 3);
    // storing for res
    std::vector<vector<double>> results;
    results.reserve(50); // change later as function input

    //loop
    for (size_t i = 1000; i <=25000; i +=500) {

        std::vector<double> price= bermud_basket_opt.PriceCall(nb_steps, i, use_antithetic, use_control_variate);
        results.push_back(std::move(price));
        cout <<  "Pricing done for " << i << " simulations" << endl;
    }

    // getting the csv
    std::string filename = "C:\\Users\\faune\\numerical-finance\\Numerical_Finance\\NumericalFinance\\Results\\bermud_scv.csv";
    WriteCSV(results, filename);
}



void EuroBasket() {
    // first we set the parameters
    int nb_assets = 3;
    double maturity = 1.0;
    double strike = 60;
    size_t nb_steps = 365;
    // size_t nb_sim = 1e4;
    vector<double> spots = {100,50,60};
    vector<double> vols = {0.10,0.25,0.16};
    double rate = 0.05;
    vector<double> weights = {0.10, 0.7, 0.2};
    vector<vector<double>> correl_mat(nb_assets, vector<double>(nb_assets, 0.1));
    for (int i = 0; i < nb_assets; i++) {
        correl_mat[i][i] = 1.0;
    }
    bool use_control_variate = false;
    bool use_antithetic = false;
    // we start the pricing
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0.,1., Unif);

    EuropeanBasketOption euro_basket_opt(nb_assets, strike, maturity, rate, spots, vols, weights, correl_mat, NormBox);
    // and we loop on the number of simuls

    // storing for res
    std::vector<vector<double>> results;
    results.reserve(50); // change later as function input

    //loop
    for (size_t i = 1000; i <=25000; i +=500) {
        std::vector<double> price= euro_basket_opt.PriceCall(nb_steps, i, use_antithetic, use_control_variate);
        results.push_back(std::move(price));
        cout <<  "Pricing done for " << i << " simulations" << endl;
    }

    // getting the csv
    std::string filename = "C:\\Users\\faune\\numerical-finance\\Numerical_Finance\\NumericalFinance\\Results\\combined_res_mc.csv";
    WriteCSV(results, filename);
}

void BasketKakutani() {
    /* Logic changes a bit here compared to other functions since we need to instantiate a different Uniform generator for
     * each given number of simulations */

    // first we set the parameters
    int nb_assets = 3;
    double maturity = 1.0;
    double strike = 60;
    size_t nb_steps = 365;
    // size_t nb_sim = 1e4;
    vector<double> spots = {100,50,60};
    vector<double> vols = {0.10,0.25,0.16};
    double rate = 0.05;
    vector<double> weights = {0.10, 0.7, 0.2};
    vector<vector<double>> correl_mat(nb_assets, vector<double>(nb_assets, 0.1));
    for (int i = 0; i < nb_assets; i++) {
        correl_mat[i][i] = 1.0;
    }
    bool use_control_variate = true;
    bool use_antithetic = false;

    // and we loop on the number of simuls

    // storing for res
    std::vector<vector<double>> results;
    results.reserve(25); // change later as function input

    //loop
    for (size_t i = 1000; i <=12500; i +=500) {
        KakutaniSequence* Kakutani = new KakutaniSequence(i, nb_assets, nb_steps);
        NormalBoxMuller* NormBox = new NormalBoxMuller(0.,1., Kakutani);
        EuropeanBasketOption european_basket_opt(nb_assets, strike, maturity, rate, spots, vols, weights, correl_mat, NormBox);

        std::vector<double> price= european_basket_opt.PriceCall(nb_steps, i, use_antithetic, use_control_variate);
        results.push_back(std::move(price));
        cout <<  "Pricing done for " << i << " simulations" << endl;
    }

    // getting the csv
    std::string filename = "C:\\Users\\mager\\Downloads\\euro.csv";
    WriteCSV(results, filename);
}
