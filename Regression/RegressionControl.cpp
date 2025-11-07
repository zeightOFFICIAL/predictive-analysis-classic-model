#include <iostream>
#include <numeric> 
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>
#include "RegressionControl.h"
#include "../Records/RecordClass.h"

namespace {
    const std::string RED = "\033[31m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string BLUE = "\033[34m";
    const std::string CYAN = "\033[36m";
    const std::string MAGENTA = "\033[35m";
    const std::string RESET = "\033[0m";
}

RegressionControl::RegressionControl(const StatisticsControl& statsControl)
    : statsControl(statsControl) {}

void RegressionControl::displayResults(const RegressionMetrics& results, 
                                       const std::string& typeLabel) const {
    std::cout << "\n" << GREEN << "=== REGRESSION RESULTS ===" << RESET << "\n";
    std::cout << typeLabel << " Price = " << std::fixed << std::setprecision(6) 
              << results.beta0 << " + " << results.beta1 << " × Gold Price\n";
    
    std::cout << "\nGoodness of Fit:\n";
    std::cout << "R-squared: " << results.R2 << " (" << (results.R2 * 100) << "%)\n";
    std::cout << "TSS: " << results.TSS << " (Total Sum of Squares)\n";
    std::cout << "ESS: " << results.ESS << " (Explained Sum of Squares)\n";
    std::cout << "RSS: " << results.RSS << " (Residual Sum of Squares)\n";
    
    double sumSqResiduals = 0.0;
    for (double r : results.residuals) {
        sumSqResiduals += r * r;
    }
    double rmse = sqrt(sumSqResiduals / results.residuals.size());
    double mae = 0.0;
    for (double r : results.residuals) {
        mae += std::abs(r);
    }
    mae /= results.residuals.size();
    
    std::cout << "RMSE: " << rmse << " (Root Mean Square Error)\n";
    std::cout << "MAE: " << mae << " (Mean Absolute Error)\n";
    
    std::cout << "\n" << YELLOW << "Economic Interpretation:" << RESET << "\n";
    std::cout << "- For every $1 increase in Gold, " << typeLabel 
              << " changes by $" << std::fixed << std::setprecision(4) << results.beta1 << "\n";
    std::cout << "- Base " << typeLabel << " price when Gold is $0: $" 
              << std::fixed << std::setprecision(2) << results.beta0 << "\n";
    std::cout << "- Model explains " << (results.R2 * 100) << "% of " 
              << typeLabel << " price variation\n";
    std::cout << "- Average prediction error: $" << std::fixed << std::setprecision(2) 
              << mae << "\n";
}

void RegressionControl::runCommodityRegression(const std::string& typeName, 
                                                const std::string& typeLabel) const {
    std::cout << CYAN << "\nRunning " << typeLabel << "-Gold Regression Analysis..." << RESET << std::endl;
    
    try {
        auto goldPrices = statsControl.getTypePrices(RecordClass::GOLD);
        auto otherPrices = statsControl.getTypePrices(typeName);
        
        if (goldPrices.empty() || otherPrices.empty()) {
            std::cout << RED << "Error: Missing price data for Gold or " << typeLabel << RESET << std::endl;
            return;
        }
        
        size_t minSize = std::min(goldPrices.size(), otherPrices.size());
        if (minSize < 10) {
            std::cout << RED << "Error: Not enough common data points! Need at least 10, got " << minSize << RESET << std::endl;
            return;
        }
        
        goldPrices.resize(minSize);
        otherPrices.resize(minSize);
        
        std::cout << "Using " << minSize << " common data points for regression\n";
        
        RegressionMetrics results = RegressionClass::calculateRegression(goldPrices, otherPrices);
        
        results.residuals.resize(minSize);
        for (size_t i = 0; i < minSize; ++i) {
            results.residuals[i] = otherPrices[i] - results.fittedValues[i];
        }
        
        displayResults(results, typeLabel);
        
        std::cout << "\nGenerate regression plot? (y/n): ";
        char plotChoice;
        std::cin >> plotChoice;
        if (tolower(plotChoice) == 'y') {
            RegressionClass::generatePlot(goldPrices, otherPrices, results, typeLabel);
            std::cout << GREEN << "✓ Regression plot generated" << RESET << std::endl;
            
            std::cout << "Generate residual plot? (y/n): ";
            std::cin >> plotChoice;
            if (tolower(plotChoice) == 'y') {
                RegressionClass::plotResiduals(goldPrices, otherPrices, results, typeLabel);
                std::cout << GREEN << "✓ Residual plot generated" << RESET << std::endl;
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Error in regression analysis: " << e.what() << RESET << "\n";
    }
}

void RegressionControl::runSilverGoldRegression() const {
    runCommodityRegression(RecordClass::SILVER, "Silver");
}

void RegressionControl::runCopperGoldRegression() const {
    runCommodityRegression(RecordClass::COPPER, "Copper");
}

void RegressionControl::runMultipleRegressionSelected(const std::vector<std::string>& selectedTypes) const {
    if (selectedTypes.empty()) {
        std::cout << RED << "Error: No commodities selected for regression!\n" << RESET;
        return;
    }
    
    if (selectedTypes.size() < 2) {
        std::cout << YELLOW << "Note: For better results, select at least 2 predictors\n" << RESET;
    }
    
    std::cout << CYAN << "\n=== MULTIPLE REGRESSION: Gold vs Selected Commodities ===" << RESET << std::endl;
    
    std::cout << "Selected commodities: ";
    for (size_t i = 0; i < selectedTypes.size(); ++i) {
        std::cout << statsControl.getTypeName(selectedTypes[i]);
        if (i < selectedTypes.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "\n\n";
    
    try {
        auto goldPrices = statsControl.getTypePrices(RecordClass::GOLD);
        if (goldPrices.empty()) {
            std::cout << RED << "Error: No gold price data available!\n" << RESET;
            return;
        }
        
        std::cout << "Gold prices data points: " << goldPrices.size() << std::endl;
        
        std::vector<std::string> predictorNames;
        std::vector<std::vector<double>> predictors;
        
        for (const auto& type : selectedTypes) {
            if (type == RecordClass::GOLD) {
                std::cout << YELLOW << "Skipping Gold (cannot use as predictor for itself)\n" << RESET;
                continue;
            }
            
            auto prices = statsControl.getTypePrices(type);
            if (prices.size() == goldPrices.size() && !prices.empty()) {
                predictorNames.push_back(statsControl.getTypeName(type));
                predictors.push_back(prices);
                std::cout << GREEN << "✓ " << std::setw(12) << std::left 
                          << statsControl.getTypeName(type) 
                          << " (" << prices.size() << " points)\n" << RESET;
            } else {
                std::cout << YELLOW << "✗ " << std::setw(12) << std::left 
                          << statsControl.getTypeName(type)
                          << " - size mismatch: " << prices.size() << " vs " << goldPrices.size() << RESET << "\n";
            }
        }
        
        if (predictors.empty()) {
            std::cout << RED << "Error: No compatible predictor data available!\n" << RESET;
            return;
        }
        
        if (predictors.size() >= goldPrices.size()) {
            std::cout << RED << "Error: Too many predictors for the number of observations!\n" << RESET;
            std::cout << "You have " << predictors.size() << " predictors but only " 
                      << goldPrices.size() << " observations.\n";
            return;
        }
        
        std::cout << CYAN << "\nFinal model: " << predictors.size() << " predictors, " 
                  << goldPrices.size() << " observations" << RESET << "\n";
        
        checkMulticollinearity(predictors, predictorNames);
        
        MultipleRegressionMetrics results = 
            RegressionClass::calculateMultipleRegression(predictors, goldPrices);
        
        displayMultipleRegressionResults(results, predictorNames);
        
        if (predictors.size() > 3) {
            std::cout << "\n" << YELLOW << "Try reduced model with only significant predictors? (y/n): " << RESET;
            char choice;
            std::cin >> choice;
            if (tolower(choice) == 'y') {
                runReducedModel(predictors, predictorNames, goldPrices, results);
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Error in multiple regression: " << e.what() << RESET << "\n";
    }
}

void RegressionControl::runReducedModel(const std::vector<std::vector<double>>& allPredictors,
                                      const std::vector<std::string>& allPredictorNames,
                                      const std::vector<double>& goldPrices,
                                      const MultipleRegressionMetrics& fullResults) const {
    
    std::vector<std::vector<double>> significantPredictors;
    std::vector<std::string> significantNames;
    
    for (size_t i = 1; i < fullResults.coefficients.size() && i - 1 < allPredictorNames.size(); ++i) {
        if (fullResults.tpValues[i] < 0.1) { 
            significantPredictors.push_back(allPredictors[i - 1]);
            significantNames.push_back(allPredictorNames[i - 1]);
        }
    }
    
    if (significantPredictors.empty()) {
        std::cout << YELLOW << "No significant predictors found for reduced model.\n" << RESET;
        return;
    }
    
    std::cout << CYAN << "\n=== REDUCED MODEL (Significant Predictors Only) ===" << RESET << std::endl;
    std::cout << "Selected predictors: ";
    for (size_t i = 0; i < significantNames.size(); ++i) {
        std::cout << significantNames[i];
        if (i < significantNames.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << "\n\n";
    
    MultipleRegressionMetrics reducedResults = 
        RegressionClass::calculateMultipleRegression(significantPredictors, goldPrices);    
    displayMultipleRegressionResults(reducedResults, significantNames);    
    compareModels(fullResults, reducedResults, allPredictorNames.size(), significantNames.size());
}

void RegressionControl::compareModels(const MultipleRegressionMetrics& fullModel,
                                    const MultipleRegressionMetrics& reducedModel,
                                    size_t fullPredictors,
                                    size_t reducedPredictors) const {
    
    std::cout << CYAN << "\n=== MODEL COMPARISON ===" << RESET << "\n";
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left << std::setw(25) << "Metric" 
              << std::setw(15) << "Full Model" 
              << std::setw(15) << "Reduced Model" 
              << std::setw(15) << "Difference" << "\n";
    
    std::cout << std::string(70, '-') << "\n";
    
    std::cout << std::left << std::setw(25) << "R-squared"
              << std::setw(15) << fullModel.R2
              << std::setw(15) << reducedModel.R2
              << std::setw(15) << (fullModel.R2 - reducedModel.R2) << "\n";
    
    std::cout << std::left << std::setw(25) << "Adjusted R-squared"
              << std::setw(15) << fullModel.adjustedR2
              << std::setw(15) << reducedModel.adjustedR2
              << std::setw(15) << (fullModel.adjustedR2 - reducedModel.adjustedR2) << "\n";
    
    std::cout << std::left << std::setw(25) << "Number of Predictors"
              << std::setw(15) << fullPredictors
              << std::setw(15) << reducedPredictors
              << std::setw(15) << (fullPredictors - reducedPredictors) << "\n";
    
    double fullRMSE = std::sqrt(fullModel.RSS / fullModel.fittedValues.size());
    double reducedRMSE = std::sqrt(reducedModel.RSS / reducedModel.fittedValues.size());
    
    std::cout << std::left << std::setw(25) << "RMSE"
              << std::setw(15) << fullRMSE
              << std::setw(15) << reducedRMSE
              << std::setw(15) << (fullRMSE - reducedRMSE) << "\n";
    
    std::cout << "\n" << YELLOW << "Interpretation:" << RESET << "\n";
    if (reducedModel.adjustedR2 > fullModel.adjustedR2) {
        std::cout << GREEN << "✓ Reduced model is better (higher adjusted R²)\n" << RESET;
    } else if (std::abs(reducedModel.adjustedR2 - fullModel.adjustedR2) < 0.01) {
        std::cout << YELLOW << "∼ Models are similar in performance\n" << RESET;
    } else {
        std::cout << "Full model has better fit but may be overfitting\n";
    }
}

void RegressionControl::runMultipleRegressionAllCommodities() const {
    std::cout << CYAN << "\n=== MULTIPLE REGRESSION: Gold vs All Commodities ===" << RESET << std::endl;
    
    try {
        auto goldPrices = statsControl.getTypePrices(RecordClass::GOLD);
        if (goldPrices.empty()) {
            std::cout << RED << "Error: No gold price data available!" << RESET << std::endl;
            return;
        }
        
        std::cout << "Gold prices data points: " << goldPrices.size() << std::endl;
        
        auto allTypes = statsControl.getAvailableTypes();
        std::vector<std::string> predictorTypes;
        std::vector<std::string> predictorNames;
        std::vector<std::vector<double>> predictors;
        
        for (const auto& type : allTypes) {
            if (type != RecordClass::GOLD) {
                auto prices = statsControl.getTypePrices(type);
                if (prices.size() == goldPrices.size() && !prices.empty()) {
                    predictorTypes.push_back(type);
                    predictorNames.push_back(statsControl.getTypeName(type));
                    predictors.push_back(prices);
                    std::cout << GREEN << "✓ " << std::setw(12) << std::left 
                              << statsControl.getTypeName(type) 
                              << " (" << prices.size() << " points)" << RESET << "\n";
                } else {
                    std::cout << YELLOW << "✗ " << std::setw(12) << std::left 
                              << statsControl.getTypeName(type)
                              << " - size mismatch: " << prices.size() << " vs " << goldPrices.size() << RESET << "\n";
                }
            }
        }
        
        if (predictors.empty()) {
            std::cout << RED << "Error: No compatible predictor data available!" << RESET << std::endl;
            return;
        }
        
        std::cout << CYAN << "\nFinal model: " << predictors.size() << " predictors, " 
                  << goldPrices.size() << " observations" << RESET << "\n";
        
        checkMulticollinearity(predictors, predictorNames);
        
        MultipleRegressionMetrics results = 
            RegressionClass::calculateMultipleRegression(predictors, goldPrices);
        
        displayMultipleRegressionResults(results, predictorNames);
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Error in multiple regression: " << e.what() << RESET << "\n";
    }
}

void RegressionControl::checkMulticollinearity(const std::vector<std::vector<double>>& predictors,
                                             const std::vector<std::string>& predictorNames) const {
    std::cout << YELLOW << "\nChecking for multicollinearity..." << RESET << "\n";
    
    size_t p = predictors.size();
    for (size_t i = 0; i < p; ++i) {
        for (size_t j = i + 1; j < p; ++j) {
            double corr = StatisticsClass::calculatePearsonCorrelation(predictors[i], predictors[j]);
            if (std::abs(corr) > 0.8) {
                std::cout << RED << "⚠ HIGH CORRELATION: " << predictorNames[i] << " vs " 
                          << predictorNames[j] << " = " << std::fixed << std::setprecision(3) 
                          << corr << RESET << "\n";
            } else if (std::abs(corr) > 0.6) {
                std::cout << YELLOW << "⚠ Moderate correlation: " << predictorNames[i] << " vs " 
                          << predictorNames[j] << " = " << std::fixed << std::setprecision(3) 
                          << corr << RESET << "\n";
            }
        }
    }
}

void RegressionControl::displayMultipleRegressionResults(const MultipleRegressionMetrics& results,
                                                       const std::vector<std::string>& predictorNames) const {
    RegressionClass::printMultipleRegressionResults(results, predictorNames);
    
    std::cout << "\n" << YELLOW << "=== ECONOMIC INTERPRETATION ===" << RESET << "\n";
    
    if (results.FpValue < 0.05) {
        std::cout << GREEN << "✓ The overall regression model is statistically significant " 
                  << "(F-test p-value: " << std::fixed << std::setprecision(4) << results.FpValue << ")" << RESET << "\n";
    } else {
        std::cout << RED << "✗ The overall regression model is NOT statistically significant " 
                  << "(F-test p-value: " << std::fixed << std::setprecision(4) << results.FpValue << ")" << RESET << "\n";
    }
    
    std::cout << "✓ The model explains " << std::fixed << std::setprecision(2) 
              << (results.R2 * 100) << "% of gold price variation\n";
    
    if (!results.coefficients.empty()) {
        std::cout << "\nIndividual coefficient interpretation:\n";
        std::cout << "- Base gold price (when all predictors are zero): $" 
                  << std::fixed << std::setprecision(2) << results.coefficients[0] << "\n";
        
        for (size_t i = 1; i < results.coefficients.size() && i - 1 < predictorNames.size(); ++i) {
            std::cout << "- " << predictorNames[i - 1] << ": ";
            std::cout << "coefficient = " << std::fixed << std::setprecision(6) << results.coefficients[i];
            
            if (results.tpValues[i] < 0.05) {
                std::cout << GREEN << " (statistically significant)" << RESET;
            } else if (results.tpValues[i] < 0.1) {
                std::cout << YELLOW << " (marginally significant)" << RESET;
            } else {
                std::cout << RED << " (not significant)" << RESET;
            }
            std::cout << "\n";
        }
    }
    
    double rmse = std::sqrt(results.RSS / results.fittedValues.size());
    double mae = 0.0;
    for (double r : results.residuals) {
        mae += std::abs(r);
    }
    mae /= results.residuals.size();
    
    std::cout << "\n" << BLUE << "Model Diagnostics:" << RESET << "\n";
    std::cout << "Root Mean Square Error (RMSE): $" << std::fixed << std::setprecision(2) << rmse << "\n";
    std::cout << "Mean Absolute Error (MAE): $" << std::fixed << std::setprecision(2) << mae << "\n";
    
    std::cout << "\n" << MAGENTA << "Significance codes: *** p < 0.001, ** p < 0.01, * p < 0.05, . p < 0.1" << RESET << "\n";
}