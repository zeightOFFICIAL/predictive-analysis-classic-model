#include "SeriesControl.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <filesystem>

SeriesControl::SeriesControl(const RecordClass& recordData, const std::string& commodityName)
    : series(createSeriesFromRecord(recordData, commodityName)) {
}

SeriesClass SeriesControl::createSeriesFromRecord(const RecordClass& recordData, const std::string& commodityName) {
    
    auto prices = recordData.getAllPrices(commodityName);
    std::vector<double> values;
    std::vector<std::string> timestamps;
    
    for (const auto& [date, price] : prices) {
        values.push_back(static_cast<double>(price));
        timestamps.push_back(recordData.formatDate(date));
    }
    
    return SeriesClass(values, timestamps, commodityName);
}

void SeriesControl::processAnomalies(double criticalValue) {
    std::cout << "=== ANOMALY PROCESSING ===" << std::endl;
    std::cout << "Original series: " << series.getName() << std::endl;
    std::cout << "Number of points: " << series.size() << std::endl;
    
    
    auto anomalies = series.detectAnomaliesIrwin(criticalValue);
    std::cout << "Anomalies detected: " << anomalies.size() << std::endl;
    
    if (!anomalies.empty()) {
        std::cout << "Anomaly indices: ";
        for (size_t idx : anomalies) {
            std::cout << idx << " ";
        }
        std::cout << std::endl;
        
        
        auto originalData = series.getData();
        
        
        series.interpolateAnomalies(anomalies);
        
        std::cout << "Anomalies successfully interpolated" << std::endl;
    }
    std::cout << std::endl;
}

SeriesClass SeriesControl::movingAverage5() const {
    auto smoothedData = series.movingAverage(5);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " (MA5)");
    return result;
}

SeriesClass SeriesControl::weightedMovingAverage5() const {
    auto weights = SeriesClass::generateLinearWeights(5);
    auto smoothedData = series.weightedMovingAverage(weights);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " (WMA5)");
    return result;
}

SeriesClass SeriesControl::weightedMovingAverage7() const {
    auto weights = SeriesClass::generateLinearWeights(7);
    auto smoothedData = series.weightedMovingAverage(weights);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " (WMA7)");
    return result;
}

SeriesClass SeriesControl::chronologicalAverage12() const {
    
    auto weights = SeriesClass::generateTriangularWeights(12);
    auto smoothedData = series.weightedMovingAverage(weights);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " (CA12)");
    return result;
}

SeriesClass SeriesControl::exponentialSmoothing(double alpha) const {
    auto smoothedData = series.exponentialSmoothing(alpha);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, 
                      series.getName() + " (Exp a=" + std::to_string(alpha).substr(0, 3) + ")");
    return result;
}

void SeriesControl::displaySeriesInfo() const {
    std::cout << "=== SERIES INFORMATION ===" << std::endl;
    std::cout << "Name: " << series.getName() << std::endl;
    std::cout << "Number of points: " << series.size() << std::endl;
    
    if (series.size() > 0) {
        auto data = series.getData();
        double min = *std::min_element(data.begin(), data.end());
        double max = *std::max_element(data.begin(), data.end());
        double sum = std::accumulate(data.begin(), data.end(), 0.0);
        double mean = sum / data.size();
        
        std::cout << "Minimum value: " << min << std::endl;
        std::cout << "Maximum value: " << max << std::endl;
        std::cout << "Mean value: " << mean << std::endl;
        std::cout << "Date range: " << series.getTimestamps().front() 
                  << " - " << series.getTimestamps().back() << std::endl;
    }
    std::cout << std::endl;
}

void SeriesControl::displayComparison(const std::vector<SeriesClass>& smoothedSeries) const {
    std::cout << "=== SERIES COMPARISON ===" << std::endl;
    std::cout << std::left << std::setw(12) << "Date" 
              << std::setw(12) << "Original";
    
    for (const auto& s : smoothedSeries) {
        std::cout << std::setw(15) << s.getName();
    }
    std::cout << std::endl;
    
    std::cout << std::string(12 + 12 + smoothedSeries.size() * 15, '-') << std::endl;
    
    auto timestamps = series.getTimestamps();
    auto originalData = series.getData();
    
    
    size_t pointsToShow = std::min(size_t(20), timestamps.size());
    
    for (size_t i = 0; i < pointsToShow; ++i) {
        std::cout << std::setw(12) << timestamps[i]
                  << std::setw(12) << std::fixed << std::setprecision(4) << originalData[i];
        
        for (const auto& s : smoothedSeries) {
            std::cout << std::setw(15) << std::fixed << std::setprecision(4) << s.getData()[i];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

std::string SeriesControl::sanitizeFilename(const std::string& filename) const {
    std::string sanitized = filename;
    std::replace(sanitized.begin(), sanitized.end(), ' ', '_');
    std::replace(sanitized.begin(), sanitized.end(), '=', '_');
    std::replace(sanitized.begin(), sanitized.end(), '(', '_');
    std::replace(sanitized.begin(), sanitized.end(), ')', '_');
    std::replace(sanitized.begin(), sanitized.end(), '/', '_');
    std::replace(sanitized.begin(), sanitized.end(), '\\', '_');
    return sanitized;
}

void SeriesControl::saveDataToFile(const SeriesClass& series, const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot create file " << filename << std::endl;
        return;
    }
    
    file << "# " << series.getName() << std::endl;
    file << "# Index Date Value" << std::endl;
    
    auto data = series.getData();
    auto timestamps = series.getTimestamps();
    
    for (size_t i = 0; i < data.size(); ++i) {
        file << i << " \"" << timestamps[i] << "\" " << std::fixed << std::setprecision(6) << data[i] << std::endl;
    }
    
    file.close();
}

void SeriesControl::cleanupDataFiles(const std::vector<std::string>& filenames) const {
    for (const auto& filename : filenames) {
        if (std::filesystem::exists(filename)) {
            std::filesystem::remove(filename);
        }
    }
}

void SeriesControl::plotAllSeries(const std::vector<SeriesClass>& smoothedSeries) const {
    std::cout << "=== GENERATING COMPARISON PLOT ===" << std::endl;
    
    
    std::filesystem::create_directories("plots/series");
    
    
    std::vector<std::string> dataFiles;
    
    std::string originalFile = "plots/series/original.dat";
    saveDataToFile(series, originalFile);
    dataFiles.push_back(originalFile);
    
    for (const auto& s : smoothedSeries) {
        std::string cleanName = sanitizeFilename(s.getName());
        std::string filename = "plots/series/" + cleanName + ".dat";
        saveDataToFile(s, filename);
        dataFiles.push_back(filename);
    }
    
    
    std::string scriptFile = "plots/series/comparison_plot.gnu";
    std::ofstream script(scriptFile);
    script << "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'" << std::endl;
    script << "set output 'plots/series/all_smoothing_methods.png'" << std::endl;
    script << "set title 'Time Series: " << series.getName() << " - Smoothing Methods Comparison'" << std::endl;
    script << "set xlabel 'Time Index'" << std::endl;
    script << "set ylabel 'Price'" << std::endl;
    script << "set grid" << std::endl;
    script << "set key outside center top horizontal" << std::endl;
    script << "set multiplot layout 3,2 title 'Smoothing Methods Comparison'" << std::endl;
    script << std::endl;
    
    
    script << "set title 'Moving Average (5 points)'" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 2 title 'Original', \\" << std::endl;
    script << "     'plots/series/" << sanitizeFilename(smoothedSeries[0].getName()) << ".dat' using 1:3 with lines lw 2 title '" 
           << smoothedSeries[0].getName() << "'" << std::endl;
    script << std::endl;
    
    
    script << "set title 'Weighted Moving Average (5 points)'" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 2 title 'Original', \\" << std::endl;
    script << "     'plots/series/" << sanitizeFilename(smoothedSeries[1].getName()) << ".dat' using 1:3 with lines lw 2 title '" 
           << smoothedSeries[1].getName() << "'" << std::endl;
    script << std::endl;
    
    
    script << "set title 'Weighted Moving Average (7 points)'" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 2 title 'Original', \\" << std::endl;
    script << "     'plots/series/" << sanitizeFilename(smoothedSeries[2].getName()) << ".dat' using 1:3 with lines lw 2 title '" 
           << smoothedSeries[2].getName() << "'" << std::endl;
    script << std::endl;
    
    
    script << "set title 'Chronological Average (12 points)'" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 2 title 'Original', \\" << std::endl;
    script << "     'plots/series/" << sanitizeFilename(smoothedSeries[3].getName()) << ".dat' using 1:3 with lines lw 2 title '" 
           << smoothedSeries[3].getName() << "'" << std::endl;
    script << std::endl;
    
    
    script << "set title 'Exponential Smoothing'" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 2 title 'Original', \\" << std::endl;
    script << "     'plots/series/" << sanitizeFilename(smoothedSeries[4].getName()) << ".dat' using 1:3 with lines lw 2 title '" 
           << smoothedSeries[4].getName() << "'" << std::endl;
    script << std::endl;
    
    
    script << "set title 'All Methods Comparison'" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 3 title 'Original', \\" << std::endl;
    for (size_t i = 0; i < smoothedSeries.size(); ++i) {
        script << "     'plots/series/" << sanitizeFilename(smoothedSeries[i].getName()) << ".dat' using 1:3 with lines lw 2 title '" 
               << smoothedSeries[i].getName() << "'";
        if (i < smoothedSeries.size() - 1) script << ", \\";
        script << std::endl;
    }
    
    script << "unset multiplot" << std::endl;
    script.close();
    
    
    std::string command = "gnuplot \"" + scriptFile + "\"";
    int result = std::system(command.c_str());
    
    if (result == 0) {
        std::cout << "Multi-panel plot saved as: plots/series/all_smoothing_methods.png" << std::endl;
    } else {
        std::cerr << "Error: GNU Plot execution failed for multi-panel plot" << std::endl;
    }
    
    
    std::string overviewScriptFile = "plots/series/overview_plot.gnu";
    std::ofstream overviewScript(overviewScriptFile);
    overviewScript << "set terminal pngcairo size 1600,900 enhanced font 'Arial,10'" << std::endl;
    overviewScript << "set output 'plots/series/smoothing_overview.png'" << std::endl;
    overviewScript << "set title 'Time Series Overview: " << series.getName() << " - All Smoothing Methods'" << std::endl;
    overviewScript << "set xlabel 'Time Index'" << std::endl;
    overviewScript << "set ylabel 'Price'" << std::endl;
    overviewScript << "set grid" << std::endl;
    overviewScript << "set key outside right top vertical box" << std::endl;
    overviewScript << "set key spacing 1.2" << std::endl;
    overviewScript << "set key width -8" << std::endl;
    overviewScript << "set key font ',9'" << std::endl;
    overviewScript << std::endl;
    
    
    overviewScript << "set style line 1 lw 3 lc rgb '#000000' dt 1" << std::endl;  
    overviewScript << "set style line 2 lw 1 lc rgb '#FF0000' dt 1" << std::endl;  
    overviewScript << "set style line 3 lw 1 lc rgb '#00FF00' dt 2" << std::endl;  
    overviewScript << "set style line 4 lw 1 lc rgb '#0000FF' dt 3" << std::endl;  
    overviewScript << "set style line 5 lw 1 lc rgb '#FF00FF' dt 4" << std::endl;  
    overviewScript << "set style line 6 lw 1 lc rgb '#FFA500' dt 5" << std::endl;  
    overviewScript << std::endl;
    
    overviewScript << "plot '" << originalFile << "' using 1:3 with lines ls 1 title 'Original Data', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[0].getName()) << ".dat' using 1:3 with lines ls 2 title '" 
                   << smoothedSeries[0].getName() << "', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[1].getName()) << ".dat' using 1:3 with lines ls 3 title '" 
                   << smoothedSeries[1].getName() << "', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[2].getName()) << ".dat' using 1:3 with lines ls 4 title '" 
                   << smoothedSeries[2].getName() << "', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[3].getName()) << ".dat' using 1:3 with lines ls 5 title '" 
                   << smoothedSeries[3].getName() << "', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[4].getName()) << ".dat' using 1:3 with lines ls 6 title '" 
                   << smoothedSeries[4].getName() << "'" << std::endl;
    
    overviewScript.close();
    
    
    std::string overviewCommand = "gnuplot \"" + overviewScriptFile + "\"";
    int overviewResult = std::system(overviewCommand.c_str());
    
    if (overviewResult == 0) {
        std::cout << "Overview plot saved as: plots/series/smoothing_overview.png" << std::endl;
    } else {
        std::cerr << "Error: GNU Plot execution failed for overview plot" << std::endl;
    }
    
    
    cleanupDataFiles(dataFiles);
    if (std::filesystem::exists(scriptFile)) {
        std::filesystem::remove(scriptFile);
    }
    if (std::filesystem::exists(overviewScriptFile)) {
        std::filesystem::remove(overviewScriptFile);
    }
}

void SeriesControl::plotIndividualComparison(const SeriesClass& smoothedSeries) const {
    std::string seriesName = smoothedSeries.getName();
    std::string cleanName = sanitizeFilename(seriesName);
    std::string filenameBase = "plots/series/individual_" + cleanName;
    
    std::string originalFile = filenameBase + "_original.dat";
    std::string smoothedFile = filenameBase + "_smoothed.dat";
    
    saveDataToFile(series, originalFile);
    saveDataToFile(smoothedSeries, smoothedFile);
    
    
    std::string scriptFile = filenameBase + ".gnu";
    std::ofstream script(scriptFile);
    script << "set terminal pngcairo size 1200,800 enhanced font 'Arial,12'" << std::endl;
    script << "set output '" << filenameBase << ".png'" << std::endl;
    script << "set title 'Time Series: " << series.getName() << " vs " << seriesName << "'" << std::endl;
    script << "set xlabel 'Time Index'" << std::endl;
    script << "set ylabel 'Price'" << std::endl;
    script << "set grid" << std::endl;
    script << "set key outside top center horizontal" << std::endl;
    script << "plot '" << originalFile << "' using 1:3 with lines lw 2 title 'Original Data', \\" << std::endl;
    script << "     '" << smoothedFile << "' using 1:3 with lines lw 2 title '" << seriesName << "'" << std::endl;
    script.close();
    
    
    std::string command = "gnuplot \"" + scriptFile + "\"";
    int result = std::system(command.c_str());
    
    if (result == 0) {
        std::cout << "Individual plot saved as: " << filenameBase << ".png" << std::endl;
    } else {
        std::cerr << "Error: GNU Plot execution failed for individual plot" << std::endl;
    }
    
    
    cleanupDataFiles({originalFile, smoothedFile});
    if (std::filesystem::exists(scriptFile)) {
        std::filesystem::remove(scriptFile);
    }
}

const SeriesClass& SeriesControl::getSeries() const {
    return series;
}