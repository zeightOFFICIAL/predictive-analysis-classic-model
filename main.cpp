#include "Models/StockPricesRecordClass.h"
#include <iostream>
#include <iomanip>
#include <ctime>

void printCommodityPrice(const StockPricesRecordClass& stockData, 
                        const std::string& commodityName,
                        const std::string& displayName,
                        time_t date) {
    float price = stockData.getPrice(commodityName, date);
    if (price >= 0) {
        std::cout << std::left << std::setw(12) << displayName 
                  << ": " << price << std::endl;
    }
}

void printDateRow(const StockPricesRecordClass& stockData, time_t date) {
    std::cout << "\nDate: " << stockData.formatDate(date) << std::endl;
    std::cout << "-------------------------" << std::endl;
    
    printCommodityPrice(stockData, StockPricesRecordClass::WTI_OIL, "WTI Oil", date);
    printCommodityPrice(stockData, StockPricesRecordClass::GOLD, "Gold", date);
    printCommodityPrice(stockData, StockPricesRecordClass::SILVER, "Silver", date);
    printCommodityPrice(stockData, StockPricesRecordClass::NATURAL_GAS, "Natural Gas", date);
    printCommodityPrice(stockData, StockPricesRecordClass::CORN, "Corn", date);
    printCommodityPrice(stockData, StockPricesRecordClass::WHEAT, "Wheat", date);
    printCommodityPrice(stockData, StockPricesRecordClass::SOYBEAN, "Soybean", date);
    printCommodityPrice(stockData, StockPricesRecordClass::COPPER, "Copper", date);
    printCommodityPrice(stockData, StockPricesRecordClass::PLATINUM, "Platinum", date);
    printCommodityPrice(stockData, StockPricesRecordClass::PALLADIUM, "Palladium", date);
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <csv_file>" << std::endl;
        return 1;
    }

    StockPricesRecordClass stockData;
    if (!stockData.loadFromCSV(argv[1])) {
        std::cerr << "Failed to load data from CSV file: " << argv[1] << std::endl;
        return 1;
    }

    std::vector<time_t> allDates = stockData.getAllDates();
    if (allDates.empty()) {
        std::cout << "No data found in the CSV file." << std::endl;
        return 0;
    }

    std::cout << "Commodity Prices Data Analysis" << std::endl;
    std::cout << "==============================" << std::endl;
    std::cout << "Total records: " << allDates.size() << std::endl;
    std::cout << "Date range: " << stockData.formatDate(allDates.front()) 
              << " to " << stockData.formatDate(allDates.back()) << std::endl;

    std::cout << "\nFirst Record:" << std::endl;
    printDateRow(stockData, allDates.front());

    std::cout << "\nLast Record:" << std::endl;
    printDateRow(stockData, allDates.back());

    return 0;
}