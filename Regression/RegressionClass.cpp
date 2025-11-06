#include "RegressionClass.h"
#include <numeric>
#include <cmath>
#include <fstream>
#include <cstdio> 
#include <stdexcept>

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

void RegressionAnalysis::plotResiduals(
    const std::vector<double>& goldPrices,
    const std::vector<double>& otherPrices,
    const RegressionMetrics& results,
    const std::string& commodityName) 
{
    if (goldPrices.size() != otherPrices.size() || goldPrices.size() != results.fittedValues.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    const std::string residualsFile = "residuals.dat";
    std::ofstream out(residualsFile);
    if (!out) {
        throw std::runtime_error("Failed to open residuals data file.");
    }

    for (size_t i = 0; i < goldPrices.size(); ++i) {
        double residual = otherPrices[i] - results.fittedValues[i];
        out << goldPrices[i] << " " << residual << "\n";
    }
    out.close();

    const std::string scriptFile = "plot_residuals.gp";
    std::ofstream script(scriptFile);
    if (!script) {
        throw std::runtime_error("Failed to create Gnuplot script.");
    }

    script << "set terminal pngcairo enhanced font 'Arial,12'\n"
           << "set output 'residuals_" << commodityName << ".png'\n"
           << "set title 'Residuals Plot (" << commodityName << " vs Gold)'\n"
           << "set xlabel 'Gold Price (USD)'\n"
           << "set ylabel 'Residuals (Actual - Predicted)'\n"
           << "set grid\n"
           << "set zeroaxis lt -1\n" 
           << "plot '" << residualsFile << "' using 1:2 with points pt 7 ps 0.5 lc rgb 'red' title 'Residuals', \\\n"
           << "     0 with lines lc rgb 'black' title 'Zero line'\n";
    script.close();

    int status = system(("gnuplot " + scriptFile).c_str());
    if (status != 0) {
        throw std::runtime_error("Gnuplot failed to execute.");
    }

    std::remove(residualsFile.c_str());
    std::remove(scriptFile.c_str());
}

void RegressionAnalysis::generatePlot(
    const std::vector<double>& goldPrices,
    const std::vector<double>& otherPrices,
    const RegressionMetrics& results,
    const std::string& commodityName) 
{
    if (goldPrices.size() != otherPrices.size() || goldPrices.size() != results.fittedValues.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    const std::string dataFilename = "regression.dat";
    std::ofstream dataFile(dataFilename);
    if (!dataFile.is_open()) {
        throw std::runtime_error("Failed to open data file for writing.");
    }

    for (size_t i = 0; i < goldPrices.size(); ++i) {
        dataFile << goldPrices[i] << " " << otherPrices[i] << " " 
                 << results.fittedValues[i] << "\n";
    }
    dataFile.close();

    const std::string scriptFilename = "plot_regression.gp";
    std::ofstream script(scriptFilename);
    if (!script.is_open()) {
        throw std::runtime_error("Failed to open Gnuplot script for writing.");
    }

    script << "set terminal pngcairo enhanced font 'Arial,12'\n"
           << "set output 'regression_" << commodityName << ".png'\n"
           << "set title '" << commodityName << " vs Gold Regression (R² = " << results.R2 << ")'\n"
           << "set xlabel 'Gold Price (USD)'\n"
           << "set ylabel '" << commodityName << " Price (USD)'\n"
           << "set grid\n"
           << "plot 'regression.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Actual', \\\n"
           << "     '' using 1:3 with lines lw 2 lc rgb 'red' title 'Fitted'\n";  // Reuse filename with ''
    script.close();

    int status = system(("gnuplot " + scriptFilename).c_str());
    if (status != 0) {
        throw std::runtime_error("Gnuplot execution failed.");
    }

    std::remove(dataFilename.c_str());
    std::remove(scriptFilename.c_str());
}