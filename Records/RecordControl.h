#ifndef RECORDCONTROL_H
#define RECORDCONTROL_H

#include "RecordClass.h"
#include <string>
#include <vector>

class RecordControl {
public:
    explicit RecordControl(const RecordClass& data);
    
    void displayAllData(bool showProgress = true) const;
    void displayDataForDate(const std::string& date) const;
    void displayDataSummary() const;
    void displayCommodityHistory(const std::string& commodityName) const;
    void displayDateRange(const std::string& startDate, const std::string& endDate) const;
    void displaySanctionsInfo() const;
    
    std::vector<std::string> getAvailableDatesFormatted() const;
    bool validateDate(const std::string& date) const;

private:
    const RecordClass& dataRef;

    void printHeader(const std::string& title) const;
    void printCommodityPrice(const std::string& commodityName, float price) const;
    void printSanctionsStatus(const std::string& date, int sanctionsValue) const;
    std::string formatPrice(float price) const;
};

#endif