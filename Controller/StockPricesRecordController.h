#ifndef STOCKPRICESRECORDCONTROLLER_H
#define STOCKPRICESRECORDCONTROLLER_H

#include "StockPricesRecordClass.h"
#include <string>
#include <vector>

class StockPricesRecordController {
public:
    explicit StockPricesRecordController(const StockPricesRecordClass& data);
    
    void displayAllData(bool showProgress = true) const;
    void displayDataForDate(const std::string& date) const;
    void displayDataSummary() const;
    void displayCommodityHistory(const std::string& commodityName) const;
    void displayDateRange(const std::string& startDate, const std::string& endDate) const;
    
    std::vector<std::string> getAvailableDatesFormatted() const;
    bool validateDate(const std::string& date) const;

private:
    const StockPricesRecordClass& dataRef;

    void printHeader(const std::string& title) const;
    void printCommodityPrice(const std::string& commodityName, float price) const;
    std::string formatPrice(float price) const;
};

#endif // STOCKPRICESRECORDCONTROLLER_H