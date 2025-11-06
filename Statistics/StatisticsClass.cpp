#include "StatisticsClass.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <map>
#include <tuple>
#include <vector>
#include <cmath>
#include <tuple>
#include <string>


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

        auto [w, p, norm] = performShapiroWilkTest(prices);
        stats.shapiroP = p;
        stats.shapiroW = w;
        stats.normality = norm;
        
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

static double normInv(double p) {
    if (p <= 0.0) return -INFINITY;
    if (p >= 1.0) return INFINITY;

    static const double a1 = -3.969683028665376e+01;
    static const double a2 =  2.209460984245205e+02;
    static const double a3 = -2.759285104469687e+02;
    static const double a4 =  1.383577518672690e+02;
    static const double a5 = -3.066479806614716e+01;
    static const double a6 =  2.506628277459239e+00;

    static const double b1 = -5.447609879822406e+01;
    static const double b2 =  1.615858368580409e+02;
    static const double b3 = -1.556989798598866e+02;
    static const double b4 =  6.680131188771972e+01;
    static const double b5 = -1.328068155288572e+01;

    static const double c1 = -7.784894002430293e-03;
    static const double c2 = -3.223964580411365e-01;
    static const double c3 = -2.400758277161838e+00;
    static const double c4 = -2.549732539343734e+00;
    static const double c5 =  4.374664141464968e+00;
    static const double c6 =  2.938163982698783e+00;

    static const double d1 =  7.784695709041462e-03;
    static const double d2 =  3.224671290700398e-01;
    static const double d3 =  2.445134137142996e+00;
    static const double d4 =  3.754408661907416e+00;

    const double p_low  = 0.02425;
    const double p_high = 1.0 - p_low;

    double q, r;
    if (p < p_low) {
        q = std::sqrt(-2.0 * std::log(p));
        return (((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
               ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
    } else if (p <= p_high) {
        q = p - 0.5;
        r = q * q;
        return (((((a1*r + a2)*r + a3)*r + a4)*r + a5)*r + a6) * q /
               (((((b1*r + b2)*r + b3)*r + b4)*r + b5)*r + 1.0);
    } else {
        q = std::sqrt(-2.0 * std::log(1.0 - p));
        return -(((((c1*q + c2)*q + c3)*q + c4)*q + c5)*q + c6) /
                ((((d1*q + d2)*q + d3)*q + d4)*q + 1.0);
    }
}

std::tuple<double, double, std::string> StatisticsClass::performShapiroWilkTest(
    const std::vector<double>& values)
{
    size_t n = values.size();
    if (n < 3) {
        return { -1.0, -1.0, "n<3 → no test" };
    }

    std::vector<double> x = values;
    std::sort(x.begin(), x.end());

    double mean = std::accumulate(x.begin(), x.end(), 0.0) / n;
    double ss = 0.0;
    for (double v : x) ss += (v - mean) * (v - mean);

    if (ss == 0.0) {
        return { 1.0, 1.0, "Perfect normal (constant)" };
    }

    std::vector<double> m(n);
    for (size_t i = 0; i < n; ++i) {
        double u = (i + 1.0 - 0.375) / (n + 0.25);
        m[i] = normInv(u);
    }

    double sum_m2 = 0.0;
    for (double v : m) sum_m2 += v * v;
    double c = 1.0 / std::sqrt(sum_m2);

    double sum_ax = 0.0;
    for (size_t i = 0; i < n; ++i) {
        sum_ax += m[i] * x[i];
    }
    double W = (sum_ax * sum_ax * c * c) / ss;
    W = std::min(1.0, std::max(0.0, W));

    double p;
    if (W > 0.98) {
        p = 0.90 + (W - 0.98) * 2.5;
    }
    else if (W > 0.95) {
        p = 0.50 + (W - 0.95) * 13.33;
    }
    else if (W > 0.90) {
        p = 0.10 + (W - 0.90) * 8.0;
    }
    else if (W > 0.85) {
        p = 0.01 + (W - 0.85) * 0.18;
    }
    else {
        p = 0.001;
    }
    
    p = std::min(1.0, std::max(0.001, p));

    std::string interp;
    if (p > 0.10) {
        interp = "Normal";
    }
    else if (p > 0.05) {
        interp = "Borderline normal";
    }
    else if (p > 0.01) {
        interp = "Moderately non-normal";
    }
    else {
        interp = "Strongly non-normal";
    }

    return { W, p, interp };
}