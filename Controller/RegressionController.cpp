#include <iostream>
#include "RegressionController.h"
#include "StockPricesRecordClass.h"

// ANSI color codes
namespace {
    const std::string RED = "\033[31m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string BLUE = "\033[34m";
    const std::string CYAN = "\033[36m";
    const std::string RESET = "\033[0m";
}

RegressionController::RegressionController(const StatisticsController& statsController)
    : statsController(statsController) {}

void RegressionController::displayResults(const RegressionMetrics& results, 
                                       const std::string& commodityLabel) const {
    std::cout << "\n" << GREEN << "Regression Results:" << RESET << "\n";
    std::cout << commodityLabel << " Price = " << results.beta0 << " + " << results.beta1 << " x Gold Price\n";
    std::cout << "TSS: " << results.TSS << " (Total Sum of Squares)\n";
    std::cout << "ESS: " << results.ESS << " (Explained Sum of Squares)\n";
    std::cout << "RSS: " << results.RSS << " (Residual Sum of Squares)\n";
    std::cout << "R-squared: " << results.R2 << "\n";
    std::cout << "\nInterpretation:\n";
    std::cout << "- For every $1 increase in Gold, " << commodityLabel << " increases by $" << results.beta1 << "\n";
    std::cout << "- Base " << commodityLabel << " price when Gold is $0: $" << results.beta0 << "\n";
    std::cout << "- Model explains " << (results.R2 * 100) << "% of " << commodityLabel << " price variation\n";
}

void RegressionController::runCommodityRegression(const std::string& commodityName, 
                                                const std::string& commodityLabel) const {
    std::cout << CYAN << "\nRunning " << commodityLabel << "-Gold Regression Analysis..." << RESET << std::endl;
    
    try {
        // Get price data
        auto goldPrices = statsController.getCommodityPrices(StockPricesRecordClass::GOLD);
        auto otherPrices = statsController.getCommodityPrices(commodityName);
        
        if (goldPrices.empty() || otherPrices.empty()) {
            std::cout << RED << "Error: Missing price data for Gold or " << commodityLabel << RESET << std::endl;
            return;
        }
        
        // Calculate regression metrics
        RegressionMetrics results = RegressionAnalysis::calculateRegression(goldPrices, otherPrices);
        displayResults(results, commodityLabel);
        
        // Ask to generate plots
        std::cout << "\nGenerate regression plot? (y/n): ";
        char plotChoice;
        std::cin >> plotChoice;
        if (tolower(plotChoice) == 'y') {
            RegressionAnalysis::generatePlot(goldPrices, otherPrices, results, commodityLabel);
            std::cout << GREEN << "\nRegression plot saved as 'regression_" << commodityLabel << ".png'" << RESET << "\n";
            
            // Ask to generate residual plot
            std::cout << "\nGenerate residual plot? (y/n): ";
            std::cin >> plotChoice;
            if (tolower(plotChoice) == 'y') {
                RegressionAnalysis::plotResiduals(goldPrices, otherPrices, results, commodityLabel);
                std::cout << GREEN << "\nResidual plot saved as 'residuals_" << commodityLabel << ".png'" << RESET << "\n";
            }
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Error in regression analysis: " << e.what() << RESET << "\n";
    }
}

void RegressionController::runSilverGoldRegression() const {
    runCommodityRegression(StockPricesRecordClass::SILVER, "Silver");
}

void RegressionController::runCopperGoldRegression() const {
    runCommodityRegression(StockPricesRecordClass::COPPER, "Copper");
}