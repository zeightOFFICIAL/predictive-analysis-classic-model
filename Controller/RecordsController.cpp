#include "RecordsController.h"
#include <iomanip>
#include <algorithm>

const std::string GREEN_COLOR = "\033[32m";
const std::string BLUE_COLOR = "\033[34m";
const std::string YELLOW_COLOR = "\033[33m";
const std::string RED_COLOR = "\033[31m";
const std::string BOLD_TEXT = "\033[1m";
const std::string RESET_TEXT = "\033[0m";

RecordsController::RecordsController(const StockPricesRecordClass& data)
    : dataRef(data) {}

void RecordsController::printHeader(const std::string& title) const {
    std::cout << "\n" << BOLD_TEXT << BLUE_COLOR << "=== " << title << " ===" << RESET_TEXT << "\n";
}

void RecordsController::printItemPrice(const std::string& name, double price) const {
    std::cout << std::left << std::setw(15) << name 
              << ": " << GREEN_COLOR << std::fixed << std::setprecision(4) 
              << price << RESET_TEXT << "\n";
}

void RecordsController::printData() const {
    printHeader("COMPLETE PRICE HISTORY DATASET");
    auto dates = dataRef.getAllDates();
    
    if (dates.empty()) {
        std::cout << YELLOW_COLOR << "No data available!" << RESET_TEXT << "\n";
        return;
    }

    for (const auto& date : dates) {
        printDataForDate(dataRef.formatDate(date));
    }
}

void RecordsController::printDataForDate(const std::string& date) const {
    std::cout << "\n" << BOLD_TEXT << "Date: " << date << RESET_TEXT << "\n";
    std::cout << std::string(date.length() + 6, '-') << "\n";
    
    printItemPrice("Oil (WIT)", 
        dataRef.getPrice(StockPricesRecordClass::OIL, date));
    printItemPrice("Gold", 
        dataRef.getPrice(StockPricesRecordClass::GOLD, date));
    printItemPrice("Silver", 
        dataRef.getPrice(StockPricesRecordClass::SILVER, date));
    printItemPrice("Natural Gas", 
        dataRef.getPrice(StockPricesRecordClass::NATURAL_GAS, date));
    printItemPrice("Corn", 
        dataRef.getPrice(StockPricesRecordClass::CORN, date));
    printItemPrice("Wheat", 
        dataRef.getPrice(StockPricesRecordClass::WHEAT, date));
    printItemPrice("Soybean", 
        dataRef.getPrice(StockPricesRecordClass::SOYBEAN, date));
    printItemPrice("Copper", 
        dataRef.getPrice(StockPricesRecordClass::COPPER, date));
    printItemPrice("Platinum", 
        dataRef.getPrice(StockPricesRecordClass::PLATINUM, date));
    printItemPrice("Palladium", 
        dataRef.getPrice(StockPricesRecordClass::PALLADIUM, date));
}

void RecordsController::printDataSummary() const {
    printHeader("DATASET SUMMARY");
    auto dates = dataRef.getAllDates();
    
    if (dates.empty()) {
        std::cout << RED_COLOR << "No data loaded!" << RESET_TEXT << "\n";
        return;
    }

    std::cout << "Records:    " << BOLD_TEXT << dates.size() << RESET_TEXT << "\n";
    std::cout << "Date Range: " << BOLD_TEXT 
              << dataRef.formatDate(dates.front()) << " to " 
              << dataRef.formatDate(dates.back()) << RESET_TEXT << "\n";
    
    auto items = dataRef.getAllCommodities();
    std::cout << "Items: " << BOLD_TEXT << items.size() << RESET_TEXT << " (";
    for (size_t i = 0; i < std::min(static_cast<size_t>(3), items.size()); ++i) {
        if (i != 0) std::cout << ", ";
        std::cout << items[i];
    }
    if (items.size() > 3) std::cout << ", ...";
    std::cout << ")\n";
}

void RecordsController::printDataForItem(const std::string& itemName) const {
    auto prices = dataRef.getAllPrices(itemName);
    if (prices.empty()) {
        std::cout << RED_COLOR << "No data available for " << itemName << RESET_TEXT << "\n";
        return;
    }

    printHeader(itemName + " PRICE HISTORY");
    std::cout << std::left << std::setw(12) << "DATE" 
              << std::setw(12) << "PRICE" << "\n";
    std::cout << std::string(24, '-') << "\n";
    
    for (const auto& [date, price] : prices) {
        std::cout << std::setw(12) << dataRef.formatDate(date)
                  << std::fixed << std::setprecision(4) 
                  << std::setw(12) << price << "\n";
    }
}