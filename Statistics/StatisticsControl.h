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
    void showCommodityAnalysis(const std::string& typeName) const;
    void showComparativeAnalysis() const;    
    std::vector<std::string> getAvailableTypes() const;
    bool typeExists(const std::string& typeName) const;
    void showGoldCorrelations() const;
    std::vector<double> getTypePrices(const std::string& commodity) const;
    std::string getTypeName(const std::string& code) const;
    void generateScatterPlotsWithGNUplot() const;
    void generateHistograms() const;
    void generateBoxplots() const;
    void generateCorrelationMatrix() const;
    
private:
    const StatisticsClass& statsRef;
    
    struct PlotData {
        std::string type;
        std::vector<std::pair<double, double>> points;
        double correlation;
    };
    
    void printSectionHeader(const std::string& title) const;
    void printStatisticRow(const std::string& label, double value, const std::string& unit = "") const;
    std::string interpretSkewness(double value) const;
    std::string interpretKurtosis(double value) const;
    std::string interpretShapiroWilk(double value) const;
    std::string formatModeOutput(const std::vector<double>& modes, bool hasSingleMode) const;
    
    std::vector<PlotData> preparePlotData() const;
    void exportPlotData(const std::vector<PlotData>& allData) const;
    void createGNUplotScript(const PlotData& data) const;    
    std::string getSanitizedTypeName(const std::string& code) const;
    void generateCorrelationMatrixPlot(const std::vector<std::vector<double>>& correlationMatrix, 
                                     const std::vector<std::string>& commodityNames) const;
};

#endif