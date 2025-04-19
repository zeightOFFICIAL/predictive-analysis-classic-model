#ifndef STATISTICSCLASS_H
#define STATISTICSCLASS_H

#include "StockPricesRecordClass.h"
#include <vector>
#include <string>
#include <map>

class StatisticsClass {
public:
    struct Statistics {
        double mean = 0.0;
        double median = 0.0;
        std::vector<double> modes;
        bool hasSingleMode = false;
        
        double variance = 0.0;
        double standardDeviation = 0.0;
        double min = 0.0;
        double max = 0.0;
        double range = 0.0;
        double iqr = 0.0;
        
        double skewness = 0.0;
        double kurtosis = 0.0;
        
        size_t count = 0;
    };

    explicit StatisticsClass(const StockPricesRecordClass& stockData);
    void calculateAll();
    
    const StockPricesRecordClass& getDataRef() const { return dataRef; }
    const Statistics& getStatistics(const std::string& commodityName) const;
    std::vector<std::string> getAvailableCommodities() const;
    
    static double calculateMean(const std::vector<double>& values);
    static double calculateMedian(std::vector<double> values);
    static double calculateVariance(const std::vector<double>& values, double mean);
    static double calculateSkewness(const std::vector<double>& values, double mean, double stdDev);
    static double calculateKurtosis(const std::vector<double>& values, double mean, double stdDev);

    std::map<std::string, double> calculateGoldCorrelations() const;
    static double calculatePearsonCorrelation(const std::vector<double>& x, const std::vector<double>& y);
private:
    const StockPricesRecordClass& dataRef;
    std::map<std::string, Statistics> statisticsMap;
    
    std::vector<double> calculateModes(const std::vector<double>& values) const;
    std::vector<double> calculateQuartiles(std::vector<double> values) const;
    double calculatePercentile(const std::vector<double>& sortedValues, double percentile) const;
};

#endif // STATISTICSCLASS_H