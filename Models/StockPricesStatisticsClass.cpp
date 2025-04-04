#include "StockPricesStatisticsClass.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>

StockPricesStatisticsClass::StockPricesStatisticsClass(const StockPricesRecordClass& stockData)
    : stockDataRef(stockData) {
    calculateAllStatistics();
}

void StockPricesStatisticsClass::calculateAllStatistics() {
    auto commodities = stockDataRef.getAllCommodities();
    
    for (const auto& commodity : commodities) {
        auto pricesMap = stockDataRef.getAllPrices(commodity);
        std::vector<double> prices;
        prices.reserve(pricesMap.size());
        
        for (const auto& [date, price] : pricesMap) {
            prices.push_back(price);
        }
        
        if (prices.empty()) continue;
        
        Statistics stats;
        stats.count = prices.size();
        
        stats.mean = calculateMean(prices);
        stats.variance = calculateVariance(prices, stats.mean);
        stats.standardDeviation = calculateStdDev(stats.variance);
        stats.median = calculateMedian(prices);
        
        auto quartiles = calculateQuartiles(prices);
        stats.iqr = quartiles[2] - quartiles[0];
        
        auto [minIt, maxIt] = std::minmax_element(prices.begin(), prices.end());
        stats.min = *minIt;
        stats.max = *maxIt;
        stats.range = stats.max - stats.min;
        
        stats.modes = calculateModes(prices);
        stats.hasSingleMode = (stats.modes.size() == 1);
        
        if (stats.standardDeviation > 0) {
            stats.skewness = calculateSkewness(prices, stats.mean, stats.standardDeviation);
            stats.kurtosis = calculateKurtosis(prices, stats.mean, stats.standardDeviation);
        } else {
            stats.skewness = 0.0;
            stats.kurtosis = -3.0;
        }
        
        statisticsMap[commodity] = stats;
    }
}

const StockPricesStatisticsClass::Statistics& StockPricesStatisticsClass::getStatistics(
    const std::string& commodityName) const {
    static const Statistics emptyStats{};
    auto it = statisticsMap.find(commodityName);
    return (it != statisticsMap.end()) ? it->second : emptyStats;
}

std::vector<std::string> StockPricesStatisticsClass::getAvailableCommodities() const {
    std::vector<std::string> commodities;
    commodities.reserve(statisticsMap.size());
    for (const auto& [commodity, _] : statisticsMap) {
        commodities.push_back(commodity);
    }
    return commodities;
}

double StockPricesStatisticsClass::calculateMean(const std::vector<double>& values) const {
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double StockPricesStatisticsClass::calculateMedian(std::vector<double> values) const {
    std::sort(values.begin(), values.end());
    const size_t n = values.size();
    return (n % 2 == 0) ? 
           (values[n/2 - 1] + values[n/2]) / 2.0 : 
           values[n/2];
}

double StockPricesStatisticsClass::calculateVariance(const std::vector<double>& values, double mean) const {
    double sum = 0.0;
    for (double val : values) {
        sum += std::pow(val - mean, 2);
    }
    return sum / values.size(); // Population variance
}

double StockPricesStatisticsClass::calculateStdDev(double variance) const {
    return std::sqrt(variance);
}

std::vector<double> StockPricesStatisticsClass::calculateQuartiles(std::vector<double> values) const {
    std::sort(values.begin(), values.end());
    return {
        calculatePercentile(values, 0.25), // Q1
        calculatePercentile(values, 0.50), // Q2 (median)
        calculatePercentile(values, 0.75)  // Q3
    };
}

std::vector<double> StockPricesStatisticsClass::calculateModes(const std::vector<double>& values) const {
    std::map<double, int> frequencyMap;
    for (double val : values) {
        frequencyMap[val]++;
    }
    
    if (frequencyMap.empty()) return {};
    
    // Find maximum frequency
    const int maxFreq = std::max_element(
        frequencyMap.begin(), frequencyMap.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; }
    )->second;
    
    // If all values are unique
    if (maxFreq == 1) return {};
    
    // Collect all modes
    std::vector<double> modes;
    for (const auto& [value, freq] : frequencyMap) {
        if (freq == maxFreq) {
            modes.push_back(value);
        }
    }
    
    return modes;
}

double StockPricesStatisticsClass::calculateSkewness(const std::vector<double>& values, 
                                                   double mean, double stdDev) const {
    double sum = 0.0;
    for (double val : values) {
        sum += std::pow((val - mean) / stdDev, 3);
    }
    return sum / values.size();
}

double StockPricesStatisticsClass::calculateKurtosis(const std::vector<double>& values,
                                                   double mean, double stdDev) const {
    double sum = 0.0;
    for (double val : values) {
        sum += std::pow((val - mean) / stdDev, 4);
    }
    return (sum / values.size()) - 3.0; // Excess kurtosis
}

double StockPricesStatisticsClass::calculatePercentile(const std::vector<double>& sortedValues, 
                                                     double percentile) const {
    if (sortedValues.empty()) return 0.0;
    if (percentile <= 0.0) return sortedValues.front();
    if (percentile >= 1.0) return sortedValues.back();
    
    const double pos = percentile * (sortedValues.size() - 1);
    const size_t lower = static_cast<size_t>(pos);
    const double frac = pos - lower;
    
    if (lower + 1 >= sortedValues.size()) {
        return sortedValues.back();
    }
    
    return sortedValues[lower] + frac * (sortedValues[lower + 1] - sortedValues[lower]);
}