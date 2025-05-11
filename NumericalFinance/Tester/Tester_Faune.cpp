//
// Created by faune on 5/9/2025.
//

#include <iostream>
#include <vector>

#include "BermudeanBasketOption.h"
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

void TestEulerPricer();
void BasicMC();
void TestEuroBasket();
void TestBermudBasket();

int main() {
//    TestEulerPricer();
//    BasicMC();
//    TestEuroBasket();
      TestBermudBasket();
}

void TestBermudBasket() {
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

    // and we loop on the number of simuls

    // storing for res
    std::vector<vector<double>> results;
    results.reserve(50); // change later as function input

    //loop
    for (size_t i = 1000; i <=25000; i +=500) {
        KakutaniSequence* Kakutani = new KakutaniSequence(i, nb_assets, nb_steps);
        NormalBoxMuller* NormBox = new NormalBoxMuller(0.,1., Kakutani);
        BermudeanBasketOption bermud_basket_opt(nb_assets, strike, maturity, rate, spots, vols, weights, correl_mat, NormBox, 3);

        std::vector<double> price= bermud_basket_opt.PriceCall(nb_steps, i, use_antithetic, use_control_variate);
        results.push_back(std::move(price));
        cout <<  "Pricing done for " << i << " simulations" << endl;
    }

    // getting the csv
    std::string filename = "C:\\Users\\faune\\numerical-finance\\Numerical_Finance\\NumericalFinance\\Results\\bermud_kakutani.csv";
    WriteCSV(results, filename);
}

void TestEuroBasket() {
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

void BasicMC() {
    // first we price the call option
    int nb_assets = 3;
    double maturity = 1.0;
    double strike = 65;
    size_t nb_steps = 365;
    size_t nb_sim = 1e4;
    vector<double> spots = {100,50,60};
    vector<double> vols = {0.10,0.25,0.16};
    double rate = 0.05;
    vector<double> weights = {0.10, 0.7, 0.2};
    vector<vector<double>> correl_mat(nb_assets, vector<double>(nb_assets, 0.1));
    for (int i = 0; i < nb_assets; i++) {
        correl_mat[i][i] = 1.0;
    }
    bool use_control_variate = false;
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0.,1., Unif);

    EuropeanBasketOption euro_basket_opt(nb_assets, strike, maturity, rate, spots, vols, weights, correl_mat, NormBox);
    std::vector<double> results= euro_basket_opt.PriceCall(nb_steps, nb_sim, false, false);

    // getting the csv
    std::string filename = "C:\\Users\\faune\\numerical-finance\\Numerical_Finance\\NumericalFinance\\Results\\one_res_mc.csv";
    // WriteCSV(results, filename);

}



void TestEulerPricer() {
    // setting MC params
    int nb_assets = 3;
    size_t nb_sim = 1e4;
    size_t nb_steps = 365;

    // inputs for pricing
    double maturity = 1.;
    double strike = 100.;
    double rate = 0.05;
    std::vector<double> spots = {95.0, 100.0, 105.0};
    std::vector<double> vols = {0.15, 0.2, 0.25};
    std::vector<double> weights = {-0.1, 0.4, 0.7};
    std::vector<std::vector<double>> correl_mat = {
        {1.0, 0.5, 0.3},
        {0.5, 1.0, 0.4},
            {0.3, 0.4, 1.0}
    };

    // starting the simulation
    std::cout << "----------Starting Monte Carlo----------" << std::endl;
    clock_t start, end;
    start = clock();

    // we then setup the random simulators
    UniformGenerator* Unif = new EcuyerCombined();
    NormalBoxMuller* NormBox = new NormalBoxMuller(0., 1., Unif); // use BrownianND instead ??

    // then we use the function coded in BlackScholesND :
    // BlackScholesND BS(NormBox, spots, strike, vols);
    // double price = BS.PriceBasketCallOption(strike, weights, maturity, correl_mat);

    // correlated BM
    // BrownianND *brownian = new BrownianND(NormBox, spots.size(), &correl_mat);
    BSEulerND BSProcess(NormBox,nb_assets ,spots, rate, vols, &correl_mat);

    //loop for MC
    double payoffs = 0.0;
    for (size_t i = 0; i < nb_sim; i++) {
        double local_payoff = 0.0;
        BSProcess.Simulate(0,maturity, nb_steps, false);
        for (size_t j=0; j< nb_assets; j++){
            local_payoff += weights[j] * BSProcess.GetPath(j)->GetValue(maturity);
        }
        payoffs += std::max<double>(local_payoff - strike, 0.0);
    }
    double price = exp(-rate*maturity) * payoffs/nb_sim;
    end = clock();
    std::cout << "The price is : " << price  << " Time : " << (end-start)*1000.0 /CLOCKS_PER_SEC << "ms" << std::endl;


}