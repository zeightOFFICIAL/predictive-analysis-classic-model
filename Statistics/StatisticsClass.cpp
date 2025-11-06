#include "StatisticsClass.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>

StatisticsClass::StatisticsClass(const RecordClass& stockData)
    : dataRef(stockData) {
    calculateAll();
}

void StatisticsClass::calculateAll() {
    auto commodities = dataRef.getAllCommodities();
    
    for (const auto& commodity : commodities) {
        auto pricesMap = dataRef.getAllPrices(commodity);
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
        stats.standardDeviation = std::sqrt(stats.variance);
        stats.median = calculateMedian(prices);
        
        auto quartiles = calculateQuartiles(prices);
        stats.iqr = quartiles[2] - quartiles[0]; // Q3 - Q1
        
        auto [minIt, maxIt] = std::minmax_element(prices.begin(), prices.end());
        stats.min = *minIt;
        stats.max = *maxIt;
        stats.range = stats.max - stats.min;
        
        stats.modes = calculateModes(prices);
        stats.hasSingleMode = (stats.modes.size() == 1);
        
        if (stats.standardDeviation > 0) {
            stats.skewness = calculateSkewness(prices, stats.mean, stats.standardDeviation);
            stats.kurtosis = calculateKurtosis(prices, stats.mean, stats.standardDeviation);
        }
        
        statisticsMap[commodity] = stats;
    }
}

const StatisticsClass::Statistics& StatisticsClass::getStatistics(
    const std::string& commodityName) const {
    static const Statistics emptyStats{};
    auto it = statisticsMap.find(commodityName);
    return (it != statisticsMap.end()) ? it->second : emptyStats;
}

std::vector<std::string> StatisticsClass::getAvailableCommodities() const {
    std::vector<std::string> commodities;
    for (const auto& [commodity, _] : statisticsMap) {
        commodities.push_back(commodity);
    }
    return commodities;
}

double StatisticsClass::calculateMean(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

double StatisticsClass::calculateMedian(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    size_t n = values.size();
    return (n % 2 == 0) ? (values[n/2 - 1] + values[n/2]) / 2.0 : values[n/2];
}

double StatisticsClass::calculateVariance(const std::vector<double>& values, double mean) {
    double sum = 0.0;
    for (double val : values) {
        sum += std::pow(val - mean, 2);
    }
    return sum / values.size(); 
}

double StatisticsClass::calculateSkewness(const std::vector<double>& values, 
                                                   double mean, double stdDev) {
    if (stdDev == 0) return 0;
    double sum = 0.0;
    for (double val : values) {
        sum += std::pow((val - mean) / stdDev, 3);
    }
    return sum / values.size();
}

double StatisticsClass::calculateKurtosis(const std::vector<double>& values,
                                                  double mean, double stdDev) {
    if (stdDev == 0) return -3; 
    double sum = 0.0;
    for (double val : values) {
        sum += std::pow((val - mean) / stdDev, 4);
    }
    return (sum / values.size()) - 3;
}

std::vector<double> StatisticsClass::calculateModes(const std::vector<double>& values) const {
    std::map<double, int> frequencyMap;
    for (double val : values) {
        frequencyMap[val]++;
    }
    
    if (frequencyMap.empty()) return {};
    
    int maxFrequency = std::max_element(
        frequencyMap.begin(), frequencyMap.end(),
        [](const auto& a, const auto& b) { return a.second < b.second; }
    )->second;
    
    if (maxFrequency == 1) return {};
    
    std::vector<double> modes;
    for (const auto& [value, freq] : frequencyMap) {
        if (freq == maxFrequency) {
            modes.push_back(value);
        }
    }
    
    return modes;
}

std::vector<double> StatisticsClass::calculateQuartiles(std::vector<double> values) const {
    std::sort(values.begin(), values.end());
    return {
        calculatePercentile(values, 0.25), // Q1
        calculatePercentile(values, 0.50), // Q2 (median)
        calculatePercentile(values, 0.75)  // Q3
    };
}

double StatisticsClass::calculatePercentile(const std::vector<double>& sortedValues, 
                                                    double percentile) const {
    if (sortedValues.empty()) return 0.0;
    if (percentile <= 0.0) return sortedValues.front();
    if (percentile >= 1.0) return sortedValues.back();
    
    double pos = percentile * (sortedValues.size() - 1);
    size_t lower = static_cast<size_t>(pos);
    double frac = pos - lower;
    
    if (lower + 1 >= sortedValues.size()) {
        return sortedValues.back();
    }
    
    return sortedValues[lower] + frac * (sortedValues[lower + 1] - sortedValues[lower]);
}

std::map<std::string, double> StatisticsClass::calculateGoldCorrelations() const {
    std::map<std::string, double> correlations;
    
    auto goldPricesMap = dataRef.getAllPrices(RecordClass::GOLD);
    std::vector<double> goldPrices;
    for (const auto& [_, price] : goldPricesMap) {
        goldPrices.push_back(price);
    }
    
    auto commodities = dataRef.getAllCommodities();
    
    for (const auto& commodity : commodities) {
        if (commodity == RecordClass::GOLD) continue; 
        
        auto pricesMap = dataRef.getAllPrices(commodity);
        std::vector<double> commodityPrices;
        
        for (const auto& [date, price] : goldPricesMap) {
            float commodityPrice = dataRef.getPrice(commodity, date);
            if (commodityPrice != -1.0f) {
                commodityPrices.push_back(commodityPrice);
            }
        }
        
        if (commodityPrices.size() == goldPrices.size() && !commodityPrices.empty()) {
            correlations[commodity] = calculatePearsonCorrelation(goldPrices, commodityPrices);
        }
    }
    
    return correlations;
}

double StatisticsClass::calculatePearsonCorrelation(
    const std::vector<double>& x, const std::vector<double>& y) {
    
    if (x.size() != y.size() || x.empty()) {
        return 0.0;
    }
    
    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
    double sum_x2 = 0.0, sum_y2 = 0.0;
    
    for (size_t i = 0; i < x.size(); ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x2 += x[i] * x[i];
        sum_y2 += y[i] * y[i];
    }
    
    double n = static_cast<double>(x.size());
    double numerator = sum_xy - (sum_x * sum_y) / n;
    double denominator_x = sqrt((sum_x2 - (sum_x * sum_x) / n));
    double denominator_y = sqrt((sum_y2 - (sum_y * sum_y) / n));
    
    if (denominator_x == 0.0 || denominator_y == 0.0) {
        return 0.0;
    }
    
    return numerator / (denominator_x * denominator_y);
}