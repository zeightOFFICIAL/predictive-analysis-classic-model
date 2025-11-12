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
    
    size_t splitPoint = data.size() / 2;
    
    double mean1 = 0.0, mean2 = 0.0;
    
    for (size_t i = 0; i < splitPoint; ++i) {
        mean1 += data[i];
    }
    mean1 /= splitPoint;
    
    for (size_t i = splitPoint; i < data.size(); ++i) {
        mean2 += data[i];
    }
    mean2 /= (data.size() - splitPoint);
    
    double var1 = 0.0, var2 = 0.0;
    
    for (size_t i = 0; i < splitPoint; ++i) {
        var1 += (data[i] - mean1) * (data[i] - mean1);
    }
    var1 /= (splitPoint - 1);
    
    for (size_t i = splitPoint; i < data.size(); ++i) {
        var2 += (data[i] - mean2) * (data[i] - mean2);
    }
    var2 /= (data.size() - splitPoint - 1);
    
    double t_statistic = (mean1 - mean2) / std::sqrt(var1/splitPoint + var2/(data.size() - splitPoint));
    
    double critical_value = 1.96;
    
    bool has_trend = std::abs(t_statistic) > critical_value;
    
    return {has_trend, t_statistic};
}

std::pair<bool, double> SeriesClass::checkTrendFosterStewart() const {
    if (data.size() < 3) {
        return {false, 0.0};
    }

    double m_i = data[0];
    double M_i = data[0]; 
    
    int u = 0;
    int d = 0; 
    int s = 0; 
    
    for (size_t i = 1; i < data.size(); ++i) {
        if (data[i] > M_i) {
            u++;
            M_i = data[i];
        } else if (data[i] < m_i) {
            d++;
            m_i = data[i];
        } else {
            s++;
        }
    }
    
    double L = u + d;
    double T = u - d;
    
    if (data.size() > 25) {
        double n = static_cast<double>(data.size());
        double mean_L = (2.0 * std::log(n) - 3.4253) / 1.5;
        double var_L = (2.0 * std::log(n) - 1.9072) / 2.25;
        
        double mean_T = 0.0;
        double var_T = (2.0 * std::log(n) - 1.9072) / 2.25;
        
        double z_L = (L - mean_L) / std::sqrt(var_L);
        double z_T = (T - mean_T) / std::sqrt(var_T);
        
        double critical_value = 1.96;
        
        bool has_trend = (std::abs(z_L) > critical_value) || (std::abs(z_T) > critical_value);
        
        return {has_trend, std::max(std::abs(z_L), std::abs(z_T))};
    } else {
        double expected_L = 2.0 * std::log(static_cast<double>(data.size()));
        bool has_trend = std::abs(L - expected_L) > expected_L * 0.5;
        
        return {has_trend, std::abs(T)};
    }
}