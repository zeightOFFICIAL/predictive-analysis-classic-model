#include "StockPricesRecordClass.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>

const std::string StockPricesRecordClass::WTI_OIL = "CL=F_closing_price";
const std::string StockPricesRecordClass::GOLD = "GC=F_closing_price";
const std::string StockPricesRecordClass::SILVER = "SI=F_closing_price";
const std::string StockPricesRecordClass::NATURAL_GAS = "NG=F_closing_price";
const std::string StockPricesRecordClass::CORN = "ZC=F_closing_price";
const std::string StockPricesRecordClass::WHEAT = "ZW=F_closing_price";
const std::string StockPricesRecordClass::SOYBEAN = "ZS=F_closing_price";
const std::string StockPricesRecordClass::COPPER = "HG=F_closing_price";
const std::string StockPricesRecordClass::PLATINUM = "PL=F_closing_price";
const std::string StockPricesRecordClass::PALLADIUM = "PA=F_closing_price";

bool StockPricesRecordClass::loadFromCSV(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) return false;

    std::string line;
    if (!std::getline(file, line)) return false;

    std::vector<std::string> commodities;
    std::istringstream headerStream(line);
    std::string token;
    
    std::getline(headerStream, token, ',');
    
    while (std::getline(headerStream, token, ',')) {
        commodities.push_back(token);
    }

    while (std::getline(file, line)) {
        std::istringstream lineStream(line);
        std::string dateStr;
        
        if (!std::getline(lineStream, dateStr, ',')) continue;
        
        time_t date = parseDate(dateStr);
        if (date == -1) continue;
        
        dates.push_back(date);
        
        for (size_t i = 0; i < commodities.size(); ++i) {
            std::string priceStr;
            if (!std::getline(lineStream, priceStr, ',')) break;
            
            try {
                float price = std::stof(priceStr);
                priceData[commodities[i]][date] = price;
            } catch (...) {
                continue;
            }
        }
    }
    
    std::sort(dates.begin(), dates.end());
    
    return true;
}

float StockPricesRecordClass::getPrice(const std::string& commodity, const std::string& date) const {
    time_t dateTime = parseDate(date);
    if (dateTime == -1) return -1.0f;
    return getPrice(commodity, dateTime);
}

float StockPricesRecordClass::getPrice(const std::string& commodity, time_t date) const {
    auto commodityIt = priceData.find(commodity);
    if (commodityIt == priceData.end()) return -1.0f;
    
    auto priceIt = commodityIt->second.find(date);
    if (priceIt == commodityIt->second.end()) return -1.0f;
    
    return priceIt->second;
}

std::map<time_t, float> StockPricesRecordClass::getAllPrices(const std::string& commodity) const {
    auto it = priceData.find(commodity);
    if (it != priceData.end()) return it->second;
    return {};
}

std::vector<std::string> StockPricesRecordClass::getAllCommodities() const {
    std::vector<std::string> commodities;
    for (const auto& pair : priceData) {
        commodities.push_back(pair.first);
    }
    return commodities;
}

std::vector<time_t> StockPricesRecordClass::getAllDates() const {
    return dates;
}

time_t StockPricesRecordClass::getEarliestDate() const {
    if (dates.empty()) return -1;
    return dates.front();
}

time_t StockPricesRecordClass::getLatestDate() const {
    if (dates.empty()) return -1;
    return dates.back();
}

std::string StockPricesRecordClass::formatDate(time_t date) const {
    std::tm* tm = localtime(&date);
    if (!tm) return "";
    
    std::ostringstream oss;
    oss << std::setw(2) << std::setfill('0') << tm->tm_mday << "-"
        << std::setw(2) << (tm->tm_mon + 1) << "-"
        << std::setw(2) << (tm->tm_year % 100);
    return oss.str();
}

time_t StockPricesRecordClass::parseDate(const std::string& dateStr) const {
    std::tm tm = {};
    char delimiter;
    std::istringstream ss(dateStr);
    
    ss >> tm.tm_mday >> delimiter;
    if (delimiter != '-') return -1;
    
    ss >> tm.tm_mon >> delimiter;
    if (delimiter != '-') return -1;
    
    ss >> tm.tm_year;
    
    // Adjust for tm structure (month 0-11, years since 1900)
    tm.tm_mon -= 1;
    if (tm.tm_year >= 0 && tm.tm_year < 100) {
        tm.tm_year += 100; // Assume 2000s
    }
    
    time_t time = mktime(&tm);
    return (time == -1) ? -1 : time;
}