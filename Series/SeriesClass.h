#ifndef SERIESCLASS_H
#define SERIESCLASS_H

#include <vector>
#include <string>

class SeriesClass {
public:
    
    struct DecompositionResult {
        std::vector<double> original;
        std::vector<double> preliminaryTrend;
        std::vector<double> primaryTrend;
        std::vector<double> secondaryTrend;
        std::vector<double> finalTrend;
        std::vector<double> preliminarySeasonal;
        std::vector<double> primarySeasonal;
        std::vector<double> secondarySeasonal;
        std::vector<double> finalSeasonal;
        std::vector<double> preliminaryResidual;
        std::vector<double> primaryResidual;
        std::vector<double> secondaryResidual;
        std::vector<double> finalResidual;
    };

private:
    std::vector<double> data;
    std::vector<std::string> timestamps;
    std::string name;

public:
    
    SeriesClass(const std::vector<double>& values, 
                const std::vector<std::string>& times, 
                const std::string& seriesName);

    
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

    
    DecompositionResult decomposeTimeSeries(int period) const;

    
    std::vector<double> smoothThreePoints() const;
    std::vector<double> calculateFirstDifferences() const;
    std::vector<double> calculateSecondDifferences() const;
    std::vector<double> calculateRelativeFirstDifferences() const;
    std::vector<double> calculateLogFirstDifferences() const;
    std::vector<double> calculateGompertzIndicator() const;
    std::vector<double> calculateLogisticIndicator() const;
    
    
    std::pair<double, double> fitLinearPolynomial() const;
    std::pair<double, double> fitExponential() const;
    std::pair<double, double> fitGompertz() const;
    std::pair<double, double> fitLogistic() const;
    
    
    std::vector<double> predictLinear(double a, double b) const;
    std::vector<double> predictExponential(double a, double b) const;   
    
    void analyzeGrowthCurveCharacteristics() const;

    std::vector<double> predictGompertz(double k, double t0) const;
    std::vector<double> predictLogistic(double k, double A) const;
};

#endif 