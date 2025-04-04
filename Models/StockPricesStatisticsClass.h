#ifndef STOCKPRICESSTATISTICSCLASS_H
#define STOCKPRICESSTATISTICSCLASS_H

#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "StockPricesRecordClass.h"

class StockPricesStatisticsClass {
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

    explicit StockPricesStatisticsClass(const StockPricesRecordClass& stockData);  
      
    void calculateAllStatistics();    
    const Statistics& getStatistics(const std::string& commodityName) const;    
    std::vector<std::string> getAvailableCommodities() const;

private:
    const StockPricesRecordClass& stockDataRef;
    std::map<std::string, Statistics> statisticsMap;
    
    double calculateMean(const std::vector<double>& values) const;
    double calculateMedian(std::vector<double> values) const;
    double calculateVariance(const std::vector<double>& values, double mean) const;
    double calculateStdDev(double variance) const;
    std::vector<double> calculateQuartiles(std::vector<double> values) const;
    std::vector<double> calculateModes(const std::vector<double>& values) const;
    double calculateSkewness(const std::vector<double>& values, double mean, double stdDev) const;
    double calculateKurtosis(const std::vector<double>& values, double mean, double stdDev) const;    
    double calculatePercentile(const std::vector<double>& sortedValues, double percentile) const;
};

#endif // STOCKPRICESSTATISTICSCLASS_H