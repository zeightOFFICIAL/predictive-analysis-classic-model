#include "StockPricesRecordController.h"
#include <iomanip>
#include <iostream>
#include <algorithm>

namespace {
    const std::string RED = "\033[31m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string BLUE = "\033[34m";
    const std::string MAGENTA = "\033[35m";
    const std::string RESET = "\033[0m";
    const std::string BOLD = "\033[1m";
}

StockPricesRecordController::StockPricesRecordController(const StockPricesRecordClass& data)
    : dataRef(data) {}

void StockPricesRecordController::printHeader(const std::string& title) const {
    std::cout << "\n" << BOLD << BLUE << "=== " << title << " ===" << RESET << "\n";
}

void StockPricesRecordController::printCommodityPrice(const std::string& commodityName, float price) const {
    std::cout << std::left << std::setw(18) << commodityName 
              << ": " << GREEN << formatPrice(price) << RESET << "\n";
}

std::string StockPricesRecordController::formatPrice(float price) const {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << price;
    return oss.str();
}

bool StockPricesRecordController::validateDate(const std::string& date) const {
    auto dates = dataRef.getAllDates();
    return std::any_of(dates.begin(), dates.end(), 
        [this, &date](time_t d) { 
            return dataRef.formatDate(d) == date; 
        });
}

void StockPricesRecordController::displayAllData(bool showProgress) const {
    auto dates = dataRef.getAllDates();
    if (dates.empty()) {
        std::cout << RED << "No data available!" << RESET << "\n";
        return;
    }

    printHeader("COMPLETE PRICE HISTORY");
    std::cout << "Total records: " << BOLD << dates.size() << RESET << "\n\n";
    
    for (size_t i = 0; i < dates.size(); ++i) {
        if (showProgress && i > 0 && i % 10 == 0) {
            std::cout << YELLOW << "Processing record " << i << " of " 
                      << dates.size() << "..." << RESET << "\n";
        }
        displayDataForDate(dataRef.formatDate(dates[i]));
    }
}

void StockPricesRecordController::displayDataForDate(const std::string& date) const {
    if (!validateDate(date)) {
        std::cout << RED << "Invalid date: " << date << RESET << "\n";
        return;
    }

    std::cout << "\n" << BOLD << "Date: " << date << RESET << "\n";
    std::cout << std::string(date.length() + 6, '-') << "\n";
    
    printCommodityPrice("WTI Oil", dataRef.getPrice(StockPricesRecordClass::WTI_OIL, date));
    printCommodityPrice("Gold", dataRef.getPrice(StockPricesRecordClass::GOLD, date));
    printCommodityPrice("Silver", dataRef.getPrice(StockPricesRecordClass::SILVER, date));
    printCommodityPrice("Natural Gas", dataRef.getPrice(StockPricesRecordClass::NATURAL_GAS, date));
    printCommodityPrice("Corn", dataRef.getPrice(StockPricesRecordClass::CORN, date));
    printCommodityPrice("Wheat", dataRef.getPrice(StockPricesRecordClass::WHEAT, date));
    printCommodityPrice("Soybean", dataRef.getPrice(StockPricesRecordClass::SOYBEAN, date));
    printCommodityPrice("Copper", dataRef.getPrice(StockPricesRecordClass::COPPER, date));
    printCommodityPrice("Platinum", dataRef.getPrice(StockPricesRecordClass::PLATINUM, date));
    printCommodityPrice("Palladium", dataRef.getPrice(StockPricesRecordClass::PALLADIUM, date));
}

void StockPricesRecordController::displayDataSummary() const {
    auto dates = dataRef.getAllDates();
    if (dates.empty()) {
        std::cout << RED << "No data loaded!" << RESET << "\n";
        return;
    }

    printHeader("DATASET SUMMARY");
    std::cout << "Records:    " << BOLD << dates.size() << RESET << "\n";
    std::cout << "Date Range: " << BOLD 
              << dataRef.formatDate(dates.front()) << " to " 
              << dataRef.formatDate(dates.back()) << RESET << "\n";
    
    auto commodities = dataRef.getAllCommodities();
    std::cout << "Commodities: " << BOLD << commodities.size() << RESET << "\n";
    
    auto getCommodityName = [](const std::string& code) {
        if (code == StockPricesRecordClass::WTI_OIL) return "Crude Oil";
        if (code == StockPricesRecordClass::GOLD) return "Gold";
        if (code == StockPricesRecordClass::SILVER) return "Silver";
        if (code == StockPricesRecordClass::NATURAL_GAS) return "Natural Gas";
        if (code == StockPricesRecordClass::CORN) return "Corn";
        if (code == StockPricesRecordClass::WHEAT) return "Wheat";
        if (code == StockPricesRecordClass::SOYBEAN) return "Soybean";
        if (code == StockPricesRecordClass::COPPER) return "Copper";
        if (code == StockPricesRecordClass::PLATINUM) return "Platinum";
        if (code == StockPricesRecordClass::PALLADIUM) return "Palladium";
        return "N/A";
    };
    
    for (const auto& commodity : commodities) {
        auto prices = dataRef.getAllPrices(commodity);
        if (!prices.empty()) {
            auto [min, max] = std::minmax_element(prices.begin(), prices.end(),
                [](const auto& a, const auto& b) { return a.second < b.second; });
            std::cout << "  - " << std::left << std::setw(12) << getCommodityName(commodity)
                      << ": " << std::fixed << std::setprecision(4)
                      << min->second << " to " << max->second << "\n";
        }
    }
}

void StockPricesRecordController::displayCommodityHistory(const std::string& commodityName) const {
    auto prices = dataRef.getAllPrices(commodityName);
    if (prices.empty()) {
        std::cout << RED << "No data available for " << commodityName << RESET << "\n";
        return;
    }

    printHeader(commodityName + " PRICE HISTORY");
    std::cout << std::left << std::setw(12) << "DATE" 
              << std::setw(12) << "PRICE" << "\n";
    std::cout << std::string(24, '-') << "\n";
    
    for (const auto& [date, price] : prices) {
        std::cout << std::setw(12) << dataRef.formatDate(date)
                  << std::setw(12) << formatPrice(price) << "\n";
    }
}

void StockPricesRecordController::displayDateRange(const std::string& startDate, const std::string& endDate) const {
    if (!validateDate(startDate) || !validateDate(endDate)) {
        std::cout << RED << "Invalid date range provided" << RESET << "\n";
        return;
    }

    printHeader("DATE RANGE: " + startDate + " to " + endDate);
    
    auto allDates = dataRef.getAllDates();
    std::vector<time_t> rangeDates;
    std::copy_if(allDates.begin(), allDates.end(), std::back_inserter(rangeDates),
        [this, &startDate, &endDate](time_t date) {
            std::string formatted = dataRef.formatDate(date);
            return formatted >= startDate && formatted <= endDate;
        });
    
    if (rangeDates.empty()) {
        std::cout << YELLOW << "No data in specified range" << RESET << "\n";
        return;
    }

    std::cout << "Found " << rangeDates.size() << " records in range\n\n";
    for (const auto& date : rangeDates) {
        displayDataForDate(dataRef.formatDate(date));
    }
}

std::vector<std::string> StockPricesRecordController::getAvailableDatesFormatted() const {
    auto dates = dataRef.getAllDates();
    std::vector<std::string> formattedDates;
    formattedDates.reserve(dates.size());
    for (const auto& date : dates) {
        formattedDates.push_back(dataRef.formatDate(date));
    }
    return formattedDates;
}