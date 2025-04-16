#include "RegressionAnalysis.h"
#include <numeric>
#include <cmath>
#include <fstream>

RegressionMetrics RegressionAnalysis::calculateRegression(const std::vector<double>& X, const std::vector<double>& Y) {
    RegressionMetrics results;
    const size_t n = X.size();
    
    const double meanX = std::accumulate(X.begin(), X.end(), 0.0) / n;
    const double meanY = std::accumulate(Y.begin(), Y.end(), 0.0) / n;
    
    double covXY = 0.0, varX = 0.0;
    for (size_t i = 0; i < n; ++i) {
        covXY += (X[i] - meanX) * (Y[i] - meanY);
        varX += pow(X[i] - meanX, 2);
    }
    results.beta1 = covXY / varX;
    results.beta0 = meanY - results.beta1 * meanX;    
    results.fittedValues.resize(n);
    results.TSS = 0.0;
    results.RSS = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        results.fittedValues[i] = results.beta0 + results.beta1 * X[i];
        double residual = Y[i] - results.fittedValues[i];
        results.TSS += pow(Y[i] - meanY, 2);
        results.RSS += pow(residual, 2);
    }
    
    results.ESS = results.TSS - results.RSS;
    results.R2 = 1.0 - (results.RSS / results.TSS);
    
    return results;
}

void RegressionAnalysis::generatePlot(const std::vector<double>& goldPrices, 
                                   const std::vector<double>& otherPrices,
                                   const RegressionMetrics& results,
                                   const std::string& commodityName) {
    std::ofstream dataFile("regression.dat");
    for (size_t i = 0; i < goldPrices.size(); ++i) {
        dataFile << goldPrices[i] << " " << otherPrices[i] << " " 
                << results.fittedValues[i] << "\n";
    }
    
    std::ofstream script("plot_regression.gp");
    script << "set terminal pngcairo enhanced font 'Arial,12'\n"
           << "set output 'regression_" << commodityName << ".png'\n"
           << "set title '" << commodityName << " vs Gold Regression (R² = " << results.R2 << ")'\n"
           << "set xlabel 'Gold Price (USD)'\n"
           << "set ylabel '" << commodityName << " Price (USD)'\n"
           << "set grid\n"
           << "plot 'regression.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Actual', \\\n"
           << "     'regression.dat' using 1:3 with lines lw 2 lc rgb 'red' title 'Fitted'\n";
    
    system("gnuplot plot_regression.gp");
}