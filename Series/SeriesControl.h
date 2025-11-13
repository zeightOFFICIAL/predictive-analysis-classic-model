#ifndef SERIESCONTROL_H
#define SERIESCONTROL_H

#include "SeriesClass.h"
#include "../Records/RecordClass.h"
#include <vector>
#include <string>

class SeriesControl {
private:
    SeriesClass series;

    // Вспомогательные методы
    std::string simplifyCommodityName(const std::string& fullName) const;
    std::string sanitizeFilename(const std::string& filename) const;
    std::string getColorByIndex(size_t index) const;
    
    void saveDataToFile(const SeriesClass& series, const std::string& filename) const;
    void saveGrowthCurveData(const std::vector<double>& values, 
                           const std::vector<std::string>& timestamps,
                           const std::string& filename, 
                           const std::string& name) const;
    void saveGrowthCurveData(const SeriesClass& series, 
                           const std::string& filename, 
                           const std::string& name) const;
    void cleanupDataFiles(const std::vector<std::string>& filenames) const;
    
    void plotGrowthCurvesComparison(const std::vector<double>& linear_pred,
                                  const std::vector<double>& exp_pred,
                                  const std::vector<double>& gomp_pred,
                                  const std::vector<double>& log_pred) const;

public:
    // Конструктор
    SeriesControl(const RecordClass& recordData, const std::string& commodityName);

    // Основные методы работы с series
    SeriesClass createSeriesFromRecord(const RecordClass& recordData, const std::string& commodityName);
    const SeriesClass& getSeries() const;
    
    // Методы обработки аномалий
    void processAnomalies(double criticalValue = 1.5);
    
    // Методы сглаживания
    SeriesClass movingAverage5() const;
    SeriesClass weightedMovingAverage5() const;
    SeriesClass weightedMovingAverage7() const;
    SeriesClass chronologicalAverage12() const;
    SeriesClass exponentialSmoothing(double alpha = 0.3) const;
    
    // Методы анализа трендов
    void analyzeTrends() const;
    
    // Методы визуализации
    void displaySeriesInfo() const;
    void displayComparison(const std::vector<SeriesClass>& smoothedSeries) const;
    void plotAllSeries(const std::vector<SeriesClass>& smoothedSeries) const;
    void plotIndividualComparison(const SeriesClass& smoothedSeries) const;
    void plotDecomposition(int period = 12) const;
    
    // Новые методы для анализа кривых роста
    void analyzeGrowthCurves() const;
};

#endif // SERIESCONTROL_H