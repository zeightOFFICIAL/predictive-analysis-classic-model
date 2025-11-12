#ifndef SERIESCLASS_H
#define SERIESCLASS_H

#include <vector>
#include <string>
#include <map>

class SeriesClass {
private:
    std::vector<double> data;
    std::vector<std::string> timestamps;
    std::string name;

public:
    SeriesClass(const std::vector<double>& values, const std::vector<std::string>& times, const std::string& seriesName);
    
    const std::vector<double>& getData() const;
    const std::vector<std::string>& getTimestamps() const;
    std::string getName() const;
    size_t size() const;
    
    void replaceData(const std::vector<double>& newData);
    std::vector<size_t> detectAnomaliesIrwin(double criticalValue = 1.5) const;
    void interpolateAnomalies(const std::vector<size_t>& anomalyIndices);
    
    std::vector<double> movingAverage(size_t window) const;
    std::vector<double> weightedMovingAverage(const std::vector<double>& weights) const;
    std::vector<double> exponentialSmoothing(double alpha) const;
    
    static std::vector<double> generateLinearWeights(size_t n);
    static std::vector<double> generateTriangularWeights(size_t n);

    std::pair<bool, double> checkTrendMeanDifferences() const;
    std::pair<bool, double> checkTrendFosterStewart() const;
};

#endif