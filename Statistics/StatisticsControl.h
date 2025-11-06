#ifndef STATISTICSCONTROL_H
#define STATISTICSCONTROL_H

#include "StatisticsClass.h"
#include <string>
#include <vector>

class StatisticsControl {
public:
    explicit StatisticsControl(const StatisticsClass& stats);
    
    void showFullReport() const;
    void showSummaryTable() const;
    void showCommodityAnalysis(const std::string& commodityName) const;
    void showComparativeAnalysis() const;
    
    std::vector<std::string> getAvailableCommodities() const;
    bool commodityExists(const std::string& commodityName) const;
    void showGoldCorrelations() const;
    void generateScatterPlotsWithGNUplot() const;
    std::vector<double> getCommodityPrices(const std::string& commodity) const;

private:
    const StatisticsClass& statsRef;

    void printSectionHeader(const std::string& title) const;
    void printStatisticRow(const std::string& label, double value, const std::string& unit = "") const;
    std::string interpretSkewness(double value) const;
    std::string interpretKurtosis(double value) const;
    std::string formatModeOutput(const std::vector<double>& modes, bool hasSingleMode) const;

    struct PlotData {
        std::string commodity;
        std::vector<std::pair<double, double>> points;
        double correlation;
    };
    
    std::vector<PlotData> preparePlotData() const;
    void exportPlotData(const std::vector<PlotData>& allData) const;
    void createGNUplotScript(const PlotData& data) const;
    std::string getCommodityName(const std::string& code) const;
    std::string getSanitizedCommodityName(const std::string& code) const;
};

#endif // STATISTICSCONTROL_H