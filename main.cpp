#include <iostream>
#include <string>
#include <iomanip>
#include <limits>

#include "Records/RecordClass.h"
#include "Records/RecordControl.h"

#include "Statistics/StatisticsClass.h"
#include "Statistics/StatisticsControl.h"

#include "Regression/RegressionControl.h"
#include "Regression/RegressionClass.h"

const std::string RED_COLOR = "\033[31m";
const std::string GREEN_COLOR = "\033[32m";
const std::string YELLOW_COLOR = "\033[33m";
const std::string BLUE_COLOR = "\033[34m";
const std::string MAGENTA_COLOR = "\033[35m";
const std::string CYAN_COLOR = "\033[36m";
const std::string RESET = "\033[0m";
const std::string BOLD = "\033[1m";

static void displayMenu() {
    std::cout << "\n\n" << GREEN_COLOR << "=== --------------------------------------------------------------------------------------- ===" << RESET;
    std::cout << "\n" << BLUE_COLOR << "====== COMMODITY ANALYSIS MENU ======" << RESET << "\n";
    std::cout << "1. Show prices for specific date\n";
    std::cout << "2. Show dataset overview\n";
    std::cout << "3. Show full statistics report\n";
    std::cout << "4. Compare commodities\n";
    
    std::cout << "\n" << GREEN_COLOR << "===== CORRELATION ANALYSIS MENU =====" << RESET << "\n";
    std::cout << "5. Analyze correlations with Gold\n";
    std::cout << "6. Generate correlation graphs\n";
    std::cout << "7. Generate correlation matrix\n";
    
    std::cout << "\n" << MAGENTA_COLOR << "======== REGRESSION ANALYSIS ========" << RESET << "\n";
    std::cout << "8. Run Silver-Gold regression\n";
    std::cout << "9. Run Copper-Gold regression\n";
    std::cout << "10. Run custom commodity regression\n";
    std::cout << "11. Run multi-regression analysis (all predictors)\n";
    std::cout << "12. Run multi-regression analysis (LNG only)\n";
    std::cout << "13. Run multi-regression analysis (plants)\n";
    std::cout << "14. Run multi-regression analysis (metals)\n";
    std::cout << "15. Run multi-regression analysis (significant predictors++)\n";
    std::cout << "16. Run multi-regression analysis (significant predictors, two most)\n";
    
    std::cout << "\n" << CYAN_COLOR << "======== DISTRIBUTION PLOTS =========" << RESET << "\n";
    std::cout << "17. Generate histograms for all commodities\n";
    std::cout << "18. Generate boxplots for all commodities\n";
    
    std::cout << "\n" << CYAN_COLOR << "========= FICTIVE PARAMETER =========" << RESET << "\n";
    std::cout << "19. Show fictive parameter details\n";
    std::cout << "20. Run multi-regression analysis (all predictors+fictive)\n";
    
    std::cout << "\n0. Exit program\n";
    std::cout << "Enter your choice (0-23): ";
}

static void displayItemsMenu(const std::vector<std::string>& items) {
    std::cout << "\n" << CYAN_COLOR << "=== AVAILABLE COMMODITIES ===" << RESET << "\n";
    for (size_t i = 0; i < items.size(); ++i) {
        std::cout << i+1 << ". " << items[i] << "\n";
    }
    std::cout << "0. Cancel\n";
    std::cout << "Select commodity: ";
}

static void runCustomRegression(const StatisticsControl& statsControl, const RegressionControl& regControl) {
    auto items = statsControl.getAvailableTypes();
    if (items.empty()) {
        std::cout << RED_COLOR << "No commodities available for analysis!\n" << RESET;
        return;
    }

    displayItemsMenu(items);
    
    int choice;
    if (!(std::cin >> choice)) {
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cout << RED_COLOR << "Invalid input! Please enter a number.\n" << RESET;
        return;
    }

    if (choice == 0) return;
    if (choice < 0 || choice > static_cast<int>(items.size())) {
        std::cout << RED_COLOR << "Invalid selection!\n" << RESET;
        return;
    }

    const std::string& selectedCommodity = items[choice-1];
    std::string commodityLabel;
    
    if (selectedCommodity == RecordClass::WTI_OIL) commodityLabel = "Crude Oil";
    else if (selectedCommodity == RecordClass::SILVER) commodityLabel = "Silver";
    else if (selectedCommodity == RecordClass::NATURAL_GAS) commodityLabel = "Natural Gas";
    else if (selectedCommodity == RecordClass::CORN) commodityLabel = "Corn";
    else if (selectedCommodity == RecordClass::WHEAT) commodityLabel = "Wheat";
    else if (selectedCommodity == RecordClass::SOYBEAN) commodityLabel = "Soybean";
    else if (selectedCommodity == RecordClass::COPPER) commodityLabel = "Copper";
    else if (selectedCommodity == RecordClass::PLATINUM) commodityLabel = "Platinum";
    else if (selectedCommodity == RecordClass::PALLADIUM) commodityLabel = "Palladium";
    else commodityLabel = selectedCommodity;

    regControl.runCommodityRegression(selectedCommodity, commodityLabel);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << RED_COLOR << "Error: Missing input file\n" << RESET;
        std::cerr << "Usage: " << argv[0] << " <data_file.csv>\n";
        return 1;
    }

    std::cout << GREEN_COLOR << "\nLoading market data from " << argv[1] << "..." << RESET;
    RecordClass data;
    if (!data.loadFromCSV(argv[1])) {
        std::cerr << RED_COLOR << "Failed to load data file!\n" << RESET;
        return 1;
    }

    RecordControl recordControl(data);
    StatisticsClass stats(data);
    StatisticsControl statsControl(stats);
    RegressionControl regControl(statsControl);

    // Define predictor sets for different regression scenarios
    const std::vector<std::string> significantPredictors_set1 = {"NG=F_closing_price"};
    const std::vector<std::string> significantPredictors_set2 = {"ZS=F_closing_price", "ZC=F_closing_price", "ZW=F_closing_price"};
    const std::vector<std::string> significantPredictors_set3 = {"SI=F_closing_price", "HG=F_closing_price", "PA=F_closing_price", "PL=F_closing_price"};
    const std::vector<std::string> significantPredictors_set4 = {"SI=F_closing_price", "NG=F_closing_price", "HG=F_closing_price", "ZS=F_closing_price"};
    const std::vector<std::string> significantPredictors_set5 = {"NG=F_closing_price", "HG=F_closing_price"};

    int choice = 0;
    std::string input;

    while (true) {
        displayMenu();
        if (!(std::cin >> choice)) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << RED_COLOR << "Invalid input! Please enter a number.\n" << RESET;
            continue;
        }
        std::cin.ignore();

        switch (choice) {
            // Basic data operations
            case 1: { 
                auto dates = recordControl.getAvailableDatesFormatted();
                if (dates.empty()) {
                    std::cout << RED_COLOR << "No market data available!\n" << RESET;
                    break;
                }
                std::cout << "Enter trading date (e.g. " << dates.front() << "): ";
                std::getline(std::cin, input);
                recordControl.displayDataForDate(input);
                break;
            }

            case 2:
                recordControl.displayDataSummary();
                break;

            case 3:
                statsControl.showFullReport();
                break;

            case 4:
                statsControl.showComparativeAnalysis();
                break;

            // Correlation analysis
            case 5:
                statsControl.showGoldCorrelations();
                break;

            case 6:
                std::cout << CYAN_COLOR << "\nGenerating scatter plots..." << RESET;
                statsControl.generateScatterPlotsWithGNUplot();
                break;

            case 7:
                std::cout << CYAN_COLOR << "\nGenerating correlation matrix..." << RESET << "\n";
                statsControl.generateCorrelationMatrix();
                break;

            // Regression analysis
            case 8:
                regControl.runSilverGoldRegression();
                break;

            case 9:
                regControl.runCopperGoldRegression();
                break;

            case 10:
                runCustomRegression(statsControl, regControl);
                break;

            case 11:
                regControl.runMultipleRegressionAllCommodities();
                break;

            case 12:
                regControl.runMultipleRegressionSelected(significantPredictors_set1);
                break;
            
            case 13:
                regControl.runMultipleRegressionSelected(significantPredictors_set2);
                break;

            case 14:
                regControl.runMultipleRegressionSelected(significantPredictors_set3);
                break;

            case 15:
                regControl.runMultipleRegressionSelected(significantPredictors_set4);
                break;

            case 16:
                regControl.runMultipleRegressionSelected(significantPredictors_set5);
                break;

            // Distribution plots
            case 17:
                std::cout << CYAN_COLOR << "\nGenerating histograms..." << RESET << "\n";
                statsControl.generateHistograms();
                break;

            case 18:
                std::cout << CYAN_COLOR << "\nGenerating boxplots..." << RESET << "\n";
                statsControl.generateBoxplots();
                break;

            // Fictive parameter
            case 19:
                recordControl.displaySanctionsInfo();
                break;

            case 20:
                regControl.runMultipleRegressionWithDummy(significantPredictors_set1, true);
                break;

            case 0:
                std::cout << GREEN_COLOR << "\nClosing market analysis...\n" << RESET;
                return 0;

            default:
                std::cout << RED_COLOR << "Invalid option! Please choose 0-23.\n" << RESET;
        }
    }
}