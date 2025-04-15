#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>  // For file operations
#include <numeric>  // For std::accumulate
#include "Models/StockPricesRecordClass.h"
#include "Controller/StockPricesRecordController.h"
#include "Models/StatisticsClass.h"
#include "Controller/StatisticsController.h"
#include <limits>
#include <math.h>
#include <vector>
#include <cmath>
#include <stdexcept>

// ANSI color codes
const std::string RED = "\033[31m";
const std::string GREEN = "\033[32m";
const std::string YELLOW = "\033[33m";
const std::string BLUE = "\033[34m";
const std::string MAGENTA = "\033[35m";
const std::string CYAN = "\033[36m";
const std::string RESET = "\033[0m";
const std::string BOLD = "\033[1m";

void displayMenu() {
    std::cout << "\n" << BLUE << "=== COMMODITY ANALYSIS MENU ===" << RESET << "\n";
    std::cout << "1. View prices for specific date\n";
    std::cout << "2. Show dataset overview\n";
    std::cout << "3. View full statistics report\n";
    std::cout << "4. Compare commodities\n";
    std::cout << "5. Analyze correlations with Gold\n";
    std::cout << "6. Generate correlation graphs\n";
    std::cout << "7. Run Silver-Gold regression\n";  // New option
    std::cout << "8. Run Copper-Gold regression\n";  // New option
    std::cout << "9. Exit program\n";
    std::cout << "Enter your choice (1-8): ";
}

struct RegressionMetrics {
    double beta0;
    double beta1;
    double TSS;  // Total Sum of Squares
    double ESS;  // Explained Sum of Squares
    double RSS;  // Residual Sum of Squares
    double R2;   // R-squared
    std::vector<double> fittedValues;
};

RegressionMetrics calculateRegression(const std::vector<double>& X, const std::vector<double>& Y) {
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

void generatePlot(const std::vector<double>& goldPrices, 
                 const std::vector<double>& silverPrices,
                 const RegressionMetrics& results) {
    std::ofstream dataFile("regression.dat");
    for (size_t i = 0; i < goldPrices.size(); ++i) {
        dataFile << goldPrices[i] << " " << silverPrices[i] << " " 
                << results.fittedValues[i] << "\n";
    }
    
    std::ofstream script("plot_regression.gp");
    script << "set terminal pngcairo enhanced font 'Arial,12'\n"
           << "set output 'regression.png'\n"
           << "set title 'Silver vs Gold Regression (R² = " << results.R2 << ")'\n"
           << "set xlabel 'Gold Price (USD)'\n"
           << "set ylabel 'x Price (USD)'\n"
           << "set grid\n"
           << "plot 'regression.dat' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Actual', \\\n"
           << "     'regression.dat' using 1:3 with lines lw 2 lc rgb 'red' title 'Fitted'\n";
    
    // Execute GNUplot
    system("gnuplot plot_regression.gp");
    std::cout << GREEN << "\nPlot saved as 'regression.png'" << RESET << "\n";
}

void runSilverGoldRegression(const StatisticsController& statsController) {
    std::cout << CYAN << "\nRunning Silver-Gold Regression Analysis..." << RESET << std::endl;
    
    try {
        // Get price data
        auto goldPrices = statsController.getCommodityPrices(StockPricesRecordClass::GOLD);
        auto silverPrices = statsController.getCommodityPrices(StockPricesRecordClass::SILVER);
        
        if (goldPrices.empty() || silverPrices.empty()) {
            std::cout << RED << "Error: Missing price data for Gold or Silver" << RESET << std::endl;
            return;
        }
        
        // Calculate regression metrics
        RegressionMetrics results = calculateRegression(goldPrices, silverPrices);
        
        // Display results
        std::cout << "\n" << GREEN << "Regression Results:" << RESET << "\n";
        std::cout << "Silver Price = " << results.beta0 << " + " << results.beta1 << " x Gold Price\n";
        std::cout << "TSS: " << results.TSS << " (Total Sum of Squares)\n";
        std::cout << "ESS: " << results.ESS << " (Explained Sum of Squares)\n";
        std::cout << "RSS: " << results.RSS << " (Residual Sum of Squares)\n";
        std::cout << "R-squared: " << results.R2 << "\n";
        std::cout << "\nInterpretation:\n";
        std::cout << "- For every $1 increase in Gold, Silver increases by $" << results.beta1 << "\n";
        std::cout << "- Base Silver price when Gold is $0: $" << results.beta0 << "\n";
        std::cout << "- Model explains " << (results.R2 * 100) << "% of Silver price variation\n";
        
        // Ask to generate plot
        std::cout << "\nGenerate regression plot? (y/n): ";
        char plotChoice;
        std::cin >> plotChoice;
        if (tolower(plotChoice) == 'y') {
            generatePlot(goldPrices, silverPrices, results);
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Error in regression analysis: " << e.what() << RESET << "\n";
    }
}

void runCopperGoldRegression(const StatisticsController& statsController) {
    std::cout << CYAN << "\nRunning Copper-Gold Regression Analysis..." << RESET << std::endl;
    
    try {
        // Get price data
        auto goldPrices = statsController.getCommodityPrices(StockPricesRecordClass::GOLD);
        auto copperPrices = statsController.getCommodityPrices(StockPricesRecordClass::COPPER);
        
        if (goldPrices.empty() || copperPrices.empty()) {
            std::cout << RED << "Error: Missing price data for Gold or Copper" << RESET << std::endl;
            return;
        }
        
        // Calculate regression metrics
        RegressionMetrics results = calculateRegression(goldPrices, copperPrices);
        
        // Display results
        std::cout << "\n" << GREEN << "Regression Results:" << RESET << "\n";
        std::cout << "Copper Price = " << results.beta0 << " + " << results.beta1 << " x Gold Price\n";
        std::cout << "TSS: " << results.TSS << " (Total Sum of Squares)\n";
        std::cout << "ESS: " << results.ESS << " (Explained Sum of Squares)\n";
        std::cout << "RSS: " << results.RSS << " (Residual Sum of Squares)\n";
        std::cout << "R-squared: " << results.R2 << "\n";
        std::cout << "\nInterpretation:\n";
        std::cout << "- For every $1 increase in Gold, Copper increases by $" << results.beta1 << "\n";
        std::cout << "- Base Copper price when Gold is $0: $" << results.beta0 << "\n";
        std::cout << "- Model explains " << (results.R2 * 100) << "% of Copper price variation\n";
        
        // Ask to generate plot
        std::cout << "\nGenerate regression plot? (y/n): ";
        char plotChoice;
        std::cin >> plotChoice;
        if (tolower(plotChoice) == 'y') {
            generatePlot(goldPrices, copperPrices, results);
        }
        
    } catch (const std::exception& e) {
        std::cerr << RED << "Error in regression analysis: " << e.what() << RESET << "\n";
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << RED << "Error: Missing input file\n" << RESET;
        std::cerr << "Usage: " << argv[0] << " <data_file.csv>\n";
        return 1;
    }

    // Load data
    std::cout << GREEN << "\nLoading market data from " << argv[1] << "..." << RESET << std::endl;
    StockPricesRecordClass data;
    if (!data.loadFromCSV(argv[1])) {
        std::cerr << RED << "Failed to load data file!\n" << RESET;
        return 1;
    }

    // Initialize controllers
    StockPricesRecordController recordController(data);
    StatisticsClass stats(data);
    StatisticsController statsController(stats);

    int choice = 0;
    std::string input;

    while (true) {
        displayMenu();
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << RED << "Invalid input! Please enter a number.\n" << RESET;
            continue;
        }
        std::cin.ignore(); // Clear newline

        switch (choice) {
            case 1: { // View specific date
                auto dates = recordController.getAvailableDatesFormatted();
                if (dates.empty()) {
                    std::cout << RED << "No market data available!\n" << RESET;
                    break;
                }
                std::cout << "Enter trading date (e.g. " << dates.front() << "): ";
                std::getline(std::cin, input);
                recordController.displayDataForDate(input);
                break;
            }

            case 2: // Dataset overview
                recordController.displayDataSummary();
                break;

            case 3: // Full statistics
                statsController.showFullReport();
                break;

            case 4: // Comparative analysis
                statsController.showComparativeAnalysis();
                break;

            case 5: // Gold correlations
                statsController.showGoldCorrelations();
                break;

            case 6: 
                std::cout << CYAN << "\nGenerating scatter plots..." << RESET << std::endl;
                statsController.generateScatterPlotsWithGNUplot();
                break;

            case 7:
                runSilverGoldRegression(statsController);
                break;

            case 8:
                runCopperGoldRegression(statsController);
            break;

            case 9: // Exit
                std::cout << GREEN << "\nClosing market analysis...\n" << RESET;
                return 0;

            default:
                std::cout << RED << "Invalid option! Please choose 1-8.\n" << RESET;
        }
    }
}