// native libs
#include <iostream>
#include <string>
#include <iomanip>
#include <limits>

// local includes
#include "Models/StockPricesRecordClass.h"
#include "Models/StatisticsClass.h"
#include "Controller/StockPricesRecordController.h"
#include "Controller/StatisticsController.h"
#include "Controller/RegressionController.h"
#include "Regression/RegressionAnalysis.h"


// console colors (ansi)
const std::string RED = "\033[31m";
const std::string GREEN = "\033[32m";
const std::string YELLOW = "\033[33m";
const std::string BLUE = "\033[34m";
const std::string MAGENTA = "\033[35m";
const std::string CYAN = "\033[36m";
const std::string RESET = "\033[0m";
const std::string BOLD = "\033[1m";

// display menu
void displayMenu() {
    std::cout << "\n" << BLUE << "=== COMMODITY ANALYSIS MENU ===" << RESET << "\n";
    std::cout << "1. View prices for specific date\n";
    std::cout << "2. Show dataset overview\n";
    std::cout << "3. View full statistics report\n";
    std::cout << "4. Compare commodities\n";
    std::cout << "5. Analyze correlations with Gold\n";
    std::cout << "6. Generate correlation graphs\n";
    std::cout << "\n" << MAGENTA << "=== REGRESSION ANALYSIS ===" << RESET << "\n";
    std::cout << "7. Run Silver-Gold regression\n";
    std::cout << "8. Run Copper-Gold regression\n";
    std::cout << "9. Run custom commodity regression\n";
    std::cout << "\n0. Exit program\n";
    std::cout << "Enter your choice (0-9): ";
}

// display commodity list
void displayCommodityMenu(const std::vector<std::string>& commodities) {
    std::cout << "\n" << CYAN << "=== AVAILABLE COMMODITIES ===" << RESET << "\n";
    for (size_t i = 0; i < commodities.size(); ++i) {
        std::cout << i+1 << ". " << commodities[i] << "\n";
    }
    std::cout << "0. Cancel\n";
    std::cout << "Select commodity: ";
}

// control function for selection of commodity
void runCustomRegression(const StatisticsController& statsController, const RegressionController& regController) {
    auto commodities = statsController.getAvailableCommodities();
    if (commodities.empty()) {
        std::cout << RED << "No commodities available for analysis!\n" << RESET;
        return;
    }

    if (commodities.empty()) {
        std::cout << YELLOW << "No commodities available for regression with Gold\n" << RESET;
        return;
    }

    displayCommodityMenu(commodities);
    
    int choice;
    if (!(std::cin >> choice)) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << RED << "Invalid input! Please enter a number.\n" << RESET;
        return;
    }

    if (choice == 0) return;
    if (choice < 0 || choice > static_cast<int>(commodities.size())) {
        std::cout << RED << "Invalid selection!\n" << RESET;
        return;
    }

    const std::string& selectedCommodity = commodities[choice-1];
    std::string commodityLabel;
    
    // Convert commodity codes to readable names
    if (selectedCommodity == StockPricesRecordClass::WTI_OIL) commodityLabel = "Crude Oil";
    else if (selectedCommodity == StockPricesRecordClass::SILVER) commodityLabel = "Silver";
    else if (selectedCommodity == StockPricesRecordClass::NATURAL_GAS) commodityLabel = "Natural Gas";
    else if (selectedCommodity == StockPricesRecordClass::CORN) commodityLabel = "Corn";
    else if (selectedCommodity == StockPricesRecordClass::WHEAT) commodityLabel = "Wheat";
    else if (selectedCommodity == StockPricesRecordClass::SOYBEAN) commodityLabel = "Soybean";
    else if (selectedCommodity == StockPricesRecordClass::COPPER) commodityLabel = "Copper";
    else if (selectedCommodity == StockPricesRecordClass::PLATINUM) commodityLabel = "Platinum";
    else if (selectedCommodity == StockPricesRecordClass::PALLADIUM) commodityLabel = "Palladium";
    else commodityLabel = selectedCommodity;

    regController.runCommodityRegression(selectedCommodity, commodityLabel);
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
    RegressionController regController(statsController);

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

            case 6: // Generate correlation graphs
                std::cout << CYAN << "\nGenerating scatter plots..." << RESET << std::endl;
                statsController.generateScatterPlotsWithGNUplot();
                break;

            case 7: // Silver-Gold regression
                regController.runSilverGoldRegression();
                break;

            case 8: // Copper-Gold regression
                regController.runCopperGoldRegression();
                break;

            case 9: // Custom commodity regression
                runCustomRegression(statsController, regController);
                break;

            case 0: // Exit
                std::cout << GREEN << "\nClosing market analysis...\n" << RESET;
                return 0;

            default:
                std::cout << RED << "Invalid option! Please choose 0-9.\n" << RESET;
        }
    }
}