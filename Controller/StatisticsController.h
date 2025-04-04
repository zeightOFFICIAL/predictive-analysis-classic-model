#ifndef STATISTICSCONTROLLER_H
#define STATISTICSCONTROLLER_H

#include "StockPricesStatisticsClass.h"
#include <string>
#include <vector>

class StatisticsController {
public:
    explicit StatisticsController(const StockPricesStatisticsClass& stats);
    
    void showFullReport() const;
    void showSummaryTable() const;
    void showCommodityAnalysis(const std::string& commodityName) const;
    void showComparativeAnalysis() const;

private:
    const StockPricesStatisticsClass& statsRef; // Reference to statistics data

    void printSectionHeader(const std::string& title) const;
    void printStatisticRow(const std::string& label, double value, const std::string& unit = "") const;
    std::string interpretSkewness(double value) const;
    std::string interpretKurtosis(double value) const;
    std::string formatModeOutput(const std::vector<double>& modes, bool hasSingleMode) const;
};

#endif // STATISTICSCONTROLLER_H