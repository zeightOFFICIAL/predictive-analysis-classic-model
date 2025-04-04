#ifndef RECORDSCONTROLLER_H
#define RECORDSCONTROLLER_H

#include "StockPricesRecordClass.h"
#include <string>
#include <vector>

class RecordsController {
public:
    explicit RecordsController(const StockPricesRecordClass& data);
    
    // Core display methods
    void displayFullDataset() const;
    void displayDateSnapshot(const std::string& date) const;
    void displayDatasetSummary() const;
    void displayCommodityHistory(const std::string& commodityName) const;
    void displayDateRangeSummary(const std::string& startDate, const std::string& endDate) const;

    // Utility methods
    std::vector<std::string> getAvailableDates() const;
    std::vector<std::string> getAvailableCommodities() const;

private:
    const StockPricesRecordClass& dataRef;

    // Formatting helpers
    void printSectionTitle(const std::string& title) const;
    void printCommodityRow(const std::string& name, double price) const;
    std::string formatPrice(double price) const;
    bool validateDate(const std::string& date) const;
};

#endif // RECORDSCONTROLLER_H