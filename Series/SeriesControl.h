#ifndef SERIESCONTROL_H
#define SERIESCONTROL_H

#include "SeriesClass.h"
#include "../Records/RecordClass.h"
#include <vector>
#include <string>

class SeriesControl {
private:
    SeriesClass series;
    
    SeriesClass createSeriesFromRecord(const RecordClass& recordData, const std::string& commodityName);
    void saveDataToFile(const SeriesClass& series, const std::string& filename) const;
    void cleanupDataFiles(const std::vector<std::string>& filenames) const;
    std::string sanitizeFilename(const std::string& filename) const; 
    
public:
    SeriesControl(const RecordClass& recordData, const std::string& commodityName);
    
    
    void processAnomalies(double criticalValue = 1.5);
    
    
    SeriesClass movingAverage5() const;
    SeriesClass weightedMovingAverage5() const;
    SeriesClass weightedMovingAverage7() const;
    SeriesClass chronologicalAverage12() const;
    SeriesClass exponentialSmoothing(double alpha = 0.3) const;
    
    
    void displaySeriesInfo() const;
    void displayComparison(const std::vector<SeriesClass>& smoothedSeries) const;
    void plotAllSeries(const std::vector<SeriesClass>& smoothedSeries) const;
    void plotIndividualComparison(const SeriesClass& smoothedSeries) const;
    
    
    const SeriesClass& getSeries() const;
};

#endif