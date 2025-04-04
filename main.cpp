#include <iostream>
#include <string>
#include "Models/StockPricesRecordClass.h"
#include "Models/StockPricesStatisticsClass.h"
#include "Controller/StockPricesController.h"
#include "Controller/StockPricesStatisticsController.h"

const std::string RED_COLOR = "\033[31m";
const std::string RESET_TEXT = "\033[0m";
const std::string GREEN_COLOR = "\033[32m";
const std::string YELLOW_COLOR = "\033[33m";

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << RED_COLOR << "Error: Missing CSV file argument\n" << RESET_TEXT;
        return 1;
    }

    std::string csvFile = argv[1];
    std::string option = (argc > 2) ? argv[2] : "--all";
    std::string dateArg = (argc > 3) ? argv[3] : "";

    try {
        StockPricesRecordClass data;
        if (!data.loadFromCSV(csvFile)) {
            throw std::runtime_error("Failed to load CSV file");
        }
        RecordsController dataPrinter(data);
        StockPricesStatisticsClass stats(data);
        StatisticsController statsPrinter(stats);

        // Process options
       
            // Complete analysis workflow
            std::cout << GREEN_COLOR << "\n=== DATA ANALYSIS STARTED ===\n" << RESET_TEXT;
            
            // Phase 1: Data Overview
            std::cout << YELLOW_COLOR << "\n[PHASE 1] DATA OVERVIEW\n" << RESET_TEXT;
            dataPrinter.printDataSummary();
            
            // Phase 2: Raw Data Sample
            std::cout << YELLOW_COLOR << "\n[PHASE 2] SAMPLE DATA\n" << RESET_TEXT;
            auto dates = data.getAllDates();
            if (!dates.empty()) {
                dataPrinter.printDataForDate(data.formatDate(dates.front()));
                if (dates.size() > 1) {
                    std::cout << "\n" << YELLOW_COLOR << "...showing first of " << dates.size() 
                              << " records...\n" << RESET_TEXT;
                }
            
            // Phase 3: Statistical Analysis
            std::cout << YELLOW_COLOR << "\n[PHASE 3] STATISTICAL ANALYSIS\n" << RESET_TEXT;
            statsPrinter.printAllStatistics();
            
            // Phase 4: Detailed Commodity Reports
            std::cout << YELLOW_COLOR << "\n[PHASE 4] DETAILED REPORTS\n" << RESET_TEXT;
            statsPrinter.printCommodityStatistics(StockPricesRecordClass::GOLD);
            statsPrinter.printCommodityStatistics(StockPricesRecordClass::OIL);
            
            std::cout << GREEN_COLOR << "\n=== ANALYSIS COMPLETED ===\n" << RESET_TEXT;
        }
        else {
            throw std::runtime_error("Invalid option");
        }
    }
    catch (const std::exception& e) {
        std::cerr << RED_COLOR << "\nError: " << e.what() << RESET_TEXT << "\n";
        displayHelp();
        return 1;
    }

    return 0;
}