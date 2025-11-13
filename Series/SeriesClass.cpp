#include "SeriesClass.h"
#include <algorithm>
#include <numeric>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <iterator>
#include <vector>

SeriesClass::SeriesClass(const std::vector<double>& values, 
                       const std::vector<std::string>& times, 
                       const std::string& seriesName)
    : data(values), timestamps(times), name(seriesName) {
    
    if (values.size() != times.size()) {
        throw std::invalid_argument("Data and timestamps must have the same size");
    }
}

const std::vector<double>& SeriesClass::getData() const {
    return data;
}

const std::vector<std::string>& SeriesClass::getTimestamps() const {
    return timestamps;
}

std::string SeriesClass::getName() const {
    return name;
}

size_t SeriesClass::size() const {
    return data.size();
}

void SeriesClass::replaceData(const std::vector<double>& newData) {
    if (newData.size() != data.size()) {
        throw std::invalid_argument("New data must have the same size as original data");
    }
    data = newData;
}

std::vector<size_t> SeriesClass::detectAnomaliesIrwin(double criticalValue) const {
    std::vector<size_t> anomalies;
    if (data.size() < 2) return anomalies;
    
    
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    
    double sumSq = 0.0;
    for (double value : data) {
        sumSq += (value - mean) * (value - mean);
    }
    double stdDev = std::sqrt(sumSq / (data.size() - 1));
    
    if (stdDev == 0) return anomalies; 
    
    std::cout << "=== IRWIN CRITERION DETAILS ===" << std::endl;
    std::cout << "Mean (μ): " << mean << std::endl;
    std::cout << "Standard Deviation (sigma): " << stdDev << std::endl;
    std::cout << "Critical Value: " << criticalValue << std::endl;
    std::cout << "--------------------------------" << std::endl;
    
    
    for (size_t i = 1; i < data.size(); ++i) {
        double difference = std::abs(data[i] - data[i-1]);
        double lambda = difference / stdDev;
        
        std::cout << "Index " << i << ": " << data[i-1] << " -> " << data[i] 
                  << " | delta = " << difference << " | lambda = " << lambda;
        
        if (lambda > criticalValue) {
            anomalies.push_back(i);
            std::cout << " <-- ANOMALY DETECTED";
        }
        std::cout << std::endl;
    }
    
    return anomalies;
}

void SeriesClass::interpolateAnomalies(const std::vector<size_t>& anomalyIndices) {
    std::cout << "=== ANOMALY INTERPOLATION ===" << std::endl;
    for (size_t idx : anomalyIndices) {
        double originalValue = data[idx];
        
        if (idx == 0) {
            
            if (data.size() > 1) {
                data[idx] = data[idx + 1];
                std::cout << "Index " << idx << ": " << originalValue << " → " << data[idx] 
                          << " (used next value)" << std::endl;
            }
        } else if (idx == data.size() - 1) {
            
            data[idx] = data[idx - 1];
            std::cout << "Index " << idx << ": " << originalValue << " → " << data[idx] 
                      << " (used previous value)" << std::endl;
        } else {
            
            data[idx] = (data[idx - 1] + data[idx + 1]) / 2.0;
            std::cout << "Index " << idx << ": " << originalValue << " > " << data[idx] 
                      << " (linear interpolation: (" << data[idx - 1] << " + " 
                      << data[idx + 1] << ") / 2)" << std::endl;
        }
    }
}

std::vector<double> SeriesClass::movingAverage(size_t window) const {
    std::vector<double> result;
    if (data.empty() || window == 0 || window > data.size()) return result;
    
    result.reserve(data.size());
    
    for (size_t i = 0; i < data.size(); ++i) {
        size_t start = (i < window/2) ? 0 : i - window/2;
        size_t end = std::min(start + window, data.size());
        
        double sum = 0.0;
        for (size_t j = start; j < end; ++j) {
            sum += data[j];
        }
        result.push_back(sum / (end - start));
    }
    
    return result;
}

std::vector<double> SeriesClass::weightedMovingAverage(const std::vector<double>& weights) const {
    std::vector<double> result;
    if (data.empty() || weights.empty()) return result;
    
    size_t window = weights.size();
    if (window > data.size()) return result;
    
    result.reserve(data.size());
    
    for (size_t i = 0; i < data.size(); ++i) {
        if (i < window/2 || i >= data.size() - window/2) {
            
            result.push_back(data[i]);
            continue;
        }
        
        double weightedSum = 0.0;
        double weightSum = 0.0;
        
        for (size_t j = 0; j < window; ++j) {
            int dataIdx = static_cast<int>(i) - static_cast<int>(window/2) + static_cast<int>(j);
            if (dataIdx >= 0 && dataIdx < static_cast<int>(data.size())) {
                weightedSum += data[dataIdx] * weights[j];
                weightSum += weights[j];
            }
        }
        
        result.push_back(weightedSum / weightSum);
    }
    
    return result;
}

std::vector<double> SeriesClass::exponentialSmoothing(double alpha) const {
    std::vector<double> result;
    if (data.empty()) return result;
    
    result.reserve(data.size());
    result.push_back(data[0]); 
    
    for (size_t i = 1; i < data.size(); ++i) {
        double smoothed = alpha * data[i] + (1 - alpha) * result[i-1];
        result.push_back(smoothed);
    }
    
    return result;
}

std::vector<double> SeriesClass::generateLinearWeights(size_t n) {
    std::vector<double> weights(n);
    double total = n * (n + 1) / 2.0; 
    
    for (size_t i = 0; i < n; ++i) {
        weights[i] = (i + 1) / total;
    }
    
    return weights;
}

std::vector<double> SeriesClass::generateTriangularWeights(size_t n) {
    std::vector<double> weights(n);
    if (n % 2 == 0) {
        
        size_t half = n / 2;
        for (size_t i = 0; i < half; ++i) {
            weights[i] = i + 1;
            weights[n - 1 - i] = i + 1;
        }
    } else {
        
        size_t center = n / 2;
        for (size_t i = 0; i <= center; ++i) {
            weights[i] = i + 1;
            weights[n - 1 - i] = i + 1;
        }
    }
    
    
    double sum = std::accumulate(weights.begin(), weights.end(), 0.0);
    for (double& weight : weights) {
        weight /= sum;
    }
    
    return weights;
}

std::pair<bool, double> SeriesClass::checkTrendMeanDifferences() const {
    if (data.size() < 3) {
        return {false, 0.0};
    }

    size_t n = data.size();
    size_t splitPoint = n / 2;
    size_t n1 = splitPoint;
    size_t n2 = n - splitPoint;

    double mean1 = 0.0, mean2 = 0.0;
    for (size_t i = 0; i < n1; ++i) mean1 += data[i];
    mean1 /= n1;
    for (size_t i = splitPoint; i < n; ++i) mean2 += data[i];
    mean2 /= n2;

    double var1 = 0.0, var2 = 0.0;
    for (size_t i = 0; i < n1; ++i) var1 += (data[i] - mean1) * (data[i] - mean1);
    var1 /= (n1 - 1);
    for (size_t i = splitPoint; i < n; ++i) var2 += (data[i] - mean2) * (data[i] - mean2);
    var2 /= (n2 - 1);

    double t_statistic = (mean1 - mean2) / std::sqrt(var1 / n1 + var2 / n2);

    double F_statistic = (var1 > var2) ? var1 / var2 : var2 / var1;
    double df1 = n1 - 1.0;
    double df2 = n2 - 1.0;

    double z = 1.96;
    double f_critical = std::exp(z * std::sqrt(2.0 * (1.0 / df1 + 1.0 / df2)));

    double t_critical = 1.96;
    bool has_trend = std::abs(t_statistic) > t_critical;
    bool variance_diff = F_statistic > f_critical;

    std::cout << "First half mean: " << mean1 << std::endl;
    std::cout << "Second half mean: " << mean2 << std::endl;
    std::cout << "First half variance: " << var1 << std::endl;
    std::cout << "Second half variance: " << var2 << std::endl;
    std::cout << "F-statistic (variance comparison): " << F_statistic << std::endl;
    std::cout << "F critical value (alpha=0.05): " << f_critical << std::endl;
    std::cout << "T-statistic: " << t_statistic << std::endl;
    std::cout << "Critical value (alpha=0.05): +/- " << t_critical << std::endl;
    std::cout << "Variance difference: " << (variance_diff ? "YES" : "NO") << std::endl;

    return {has_trend, t_statistic};
}

std::pair<bool, double> SeriesClass::checkTrendFosterStewart() const {
    if (data.size() < 3) {
        return {false, 0.0};
    }

    size_t n = data.size();
    double m_i = data[0];
    double M_i = data[0];
    int u = 0, d = 0;

    for (size_t i = 1; i < n; ++i) {
        if (data[i] > M_i) {
            ++u;
            M_i = data[i];
        } else if (data[i] < m_i) {
            ++d;
            m_i = data[i];
        }
    }

    double S = static_cast<double>(u + d);
    double D = static_cast<double>(u - d);

    double E_S = 0.0;
    double Var_S = 0.0;
    for (size_t i = 1; i <= n; ++i) {
        E_S += 2.0 / i;
        Var_S += 4.0 / i - 8.0 / (i * i);
    }
    E_S -= 1.0;
    Var_S = std::max(Var_S, 1e-9);

    double ts = (S - E_S) / std::sqrt(Var_S);
    double td = D / std::sqrt(Var_S);

    double critical_value = 1.96;
    bool has_trend = (std::abs(ts) > critical_value) || (std::abs(td) > critical_value);

    std::cout << "S (trend score): " << S << std::endl;
    std::cout << "D (difference score): " << D << std::endl;
    std::cout << "ts (t-statistic for S): " << ts << std::endl;
    std::cout << "td (t-statistic for D): " << td << std::endl;
    std::cout << "New maxima (u): " << u << std::endl;
    std::cout << "New minima (d): " << d << std::endl;
    std::cout << "Critical value (alpha=0.05): +/- " << critical_value << std::endl;

    return {has_trend, std::max(std::abs(ts), std::abs(td))};
}

SeriesClass::DecompositionResult SeriesClass::decomposeTimeSeries(int period) const {
    DecompositionResult result;
    size_t n = data.size();
    
    if (n < period * 2) {
        throw std::invalid_argument("Series too short for decomposition");
    }
    
    result.original = data;
    result.preliminaryTrend.resize(n, 0.0);
    result.primaryTrend.resize(n, 0.0);
    result.secondaryTrend.resize(n, 0.0);
    result.finalTrend.resize(n, 0.0);
    result.preliminarySeasonal.resize(n, 0.0);
    result.primarySeasonal.resize(n, 0.0);
    result.secondarySeasonal.resize(n, 0.0);
    result.finalSeasonal.resize(n, 0.0);
    result.preliminaryResidual.resize(n, 0.0);
    result.primaryResidual.resize(n, 0.0);
    result.secondaryResidual.resize(n, 0.0);
    result.finalResidual.resize(n, 0.0);
        
    std::cout << "Decomposition info: " << n << " points, period = " << period << std::endl;
    
    int window1 = 30;
    for (size_t i = 0; i < n; ++i) {
        size_t start = (i < window1/2) ? 0 : i - window1/2;
        size_t end = std::min(start + window1, n);
        double sum = 0.0;
        size_t count = 0;
        for (size_t j = start; j < end; ++j) {
            sum += data[j];
            count++;
        }
        result.preliminaryTrend[i] = (count > 0) ? sum / count : data[i];
    }
    
    auto weights7 = generateLinearWeights(7);
    auto wma7Data = weightedMovingAverage(weights7);
    if (wma7Data.size() == n) {
        result.primaryTrend = wma7Data;
    } else {
        result.primaryTrend = result.preliminaryTrend;
    }
    
    int trendWindow = period;
    for (size_t i = 0; i < n; ++i) {
        int halfWindow = trendWindow / 2;
        int start = static_cast<int>(i) - halfWindow;
        int end = static_cast<int>(i) + halfWindow + (trendWindow % 2);
        
        if (start < 0) start = 0;
        if (end > static_cast<int>(n)) end = n;
        
        double sum = 0.0;
        int count = 0;
        for (int j = start; j < end; ++j) {
            sum += data[j];
            count++;
        }
        result.finalTrend[i] = (count > 0) ? sum / count : data[i];
    }
    
    result.secondaryTrend = result.finalTrend;
    
    std::vector<double> detrended(n, 0.0);
    for (size_t i = 0; i < n; ++i) {
        detrended[i] = data[i] - result.finalTrend[i];
    }
    
    std::vector<std::vector<double>> seasonalGroups(period);
    for (size_t i = 0; i < n; ++i) {
        int seasonIndex = i % period;
        seasonalGroups[seasonIndex].push_back(detrended[i]);
    }
    
    std::vector<double> seasonalPattern(period, 0.0);
    for (int i = 0; i < period; ++i) {
        if (!seasonalGroups[i].empty()) {
            std::vector<double> sorted = seasonalGroups[i];
            std::sort(sorted.begin(), sorted.end());
            size_t mid = sorted.size() / 2;
            if (sorted.size() % 2 == 0) {
                seasonalPattern[i] = (sorted[mid-1] + sorted[mid]) / 2.0;
            } else {
                seasonalPattern[i] = sorted[mid];
            }
        }
    }
    
    double seasonalSum = std::accumulate(seasonalPattern.begin(), seasonalPattern.end(), 0.0);
    double adjustment = seasonalSum / period;
    for (int i = 0; i < period; ++i) {
        seasonalPattern[i] -= adjustment;
    }
    
    for (size_t i = 0; i < n; ++i) {
        int seasonIndex = i % period;
        result.finalSeasonal[i] = seasonalPattern[seasonIndex];
        result.secondarySeasonal[i] = result.finalSeasonal[i];
    }
        
    for (size_t i = 0; i < n; ++i) {
        result.finalResidual[i] = data[i] - result.finalTrend[i] - result.finalSeasonal[i];
        result.secondaryResidual[i] = result.finalResidual[i];
        
        result.preliminarySeasonal[i] = data[i] - result.preliminaryTrend[i];
        result.primarySeasonal[i] = data[i] - result.primaryTrend[i];
        
        result.preliminaryResidual[i] = data[i] - (result.preliminaryTrend[i] + result.preliminarySeasonal[i]);
        result.primaryResidual[i] = data[i] - (result.primaryTrend[i] + result.primarySeasonal[i]);
    }
    
    std::cout << "Seasonal pattern statistics:" << std::endl;
    auto minMax = std::minmax_element(seasonalPattern.begin(), seasonalPattern.end());
    std::cout << "Range: [" << *minMax.first << ", " << *minMax.second << "]" << std::endl;
    std::cout << "Amplitude: " << (*minMax.second - *minMax.first) << std::endl;
    
    return result;
}