#include "SeriesControl.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <filesystem>
#include <math.h>

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
    
    std::string simplifiedName = simplifyCommodityName(commodityName);
    
    return SeriesClass(values, timestamps, simplifiedName);
}

std::string SeriesControl::simplifyCommodityName(const std::string& fullName) const {
    if (fullName == RecordClass::WTI_OIL) return "Crude Oil";
    else if (fullName == RecordClass::GOLD) return "Gold";
    else if (fullName == RecordClass::SILVER) return "Silver";
    else if (fullName == RecordClass::NATURAL_GAS) return "Natural Gas";
    else if (fullName == RecordClass::CORN) return "Corn";
    else if (fullName == RecordClass::WHEAT) return "Wheat";
    else if (fullName == RecordClass::SOYBEAN) return "Soybean";
    else if (fullName == RecordClass::COPPER) return "Copper";
    else if (fullName == RecordClass::PLATINUM) return "Platinum";
    else if (fullName == RecordClass::PALLADIUM) return "Palladium";
    else return fullName;
}

void SeriesControl::processAnomalies(double criticalValue) {
    std::cout << "=== ANOMALY PROCESSING ===" << std::endl;
    std::cout << "Original series: " << series.getName() << std::endl;
    std::cout << "Number of points: " << series.size() << std::endl;
    
    auto anomalies = series.detectAnomaliesIrwin(criticalValue);
    
    auto data = series.getData();
    if (data.size() >= 2) {
        double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
        double sumSq = 0.0;
        for (double value : data) {
            sumSq += (value - mean) * (value - mean);
        }
        double sigma = std::sqrt(sumSq / (data.size() - 1));
        
        std::cout << "--------------------------------" << std::endl;
        std::cout << "FINAL SUMMARY:" << std::endl;
        std::cout << "Standard Deviation (sigma): " << sigma << std::endl;
        std::cout << "Critical Value: " << criticalValue << std::endl;
        std::cout << "Anomalies detected: " << anomalies.size() << std::endl;
        std::cout << "=================================" << std::endl;
    }
    
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
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " - MA5");
    return result;
}

SeriesClass SeriesControl::weightedMovingAverage5() const {
    auto weights = SeriesClass::generateLinearWeights(5);
    auto smoothedData = series.weightedMovingAverage(weights);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " - WMA5");
    return result;
}

SeriesClass SeriesControl::weightedMovingAverage7() const {
    auto weights = SeriesClass::generateLinearWeights(7);
    auto smoothedData = series.weightedMovingAverage(weights);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " - WMA7");
    return result;
}

SeriesClass SeriesControl::chronologicalAverage12() const {
    auto weights = SeriesClass::generateTriangularWeights(12);
    auto smoothedData = series.weightedMovingAverage(weights);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " - CA12");
    return result;
}

SeriesClass SeriesControl::exponentialSmoothing(double alpha) const {
    auto smoothedData = series.exponentialSmoothing(alpha);
    auto timestamps = series.getTimestamps();
    
    SeriesClass result(smoothedData, timestamps, series.getName() + " - Exp");
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
    
    std::vector<std::string> residualFiles;
    auto originalData = series.getData();
    
    for (size_t i = 0; i < smoothedSeries.size(); ++i) {
        const auto& smoothed = smoothedSeries[i];
        auto smoothedData = smoothed.getData();
        std::string residualFilename = "plots/series/residual_" + sanitizeFilename(smoothed.getName()) + ".dat";
        
        std::ofstream residualFile(residualFilename);
        if (residualFile.is_open()) {
            residualFile << "# Residuals: " << smoothed.getName() << std::endl;
            residualFile << "# Index Date Residual" << std::endl;
            
            auto timestamps = series.getTimestamps();
            for (size_t j = 0; j < originalData.size() && j < smoothedData.size(); ++j) {
                double residual = originalData[j] - smoothedData[j];
                residualFile << j << " \"" << timestamps[j] << "\" " << std::fixed << std::setprecision(6) << residual << std::endl;
            }
            residualFile.close();
            residualFiles.push_back(residualFilename);
        }
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
    overviewScript << "set terminal pngcairo size 6400,3600 enhanced font 'Arial,14'" << std::endl;
    overviewScript << "set output 'plots/series/smoothing_overview.png'" << std::endl;
    overviewScript << "set multiplot layout 2,1 title 'Time Series Overview: " << series.getName() << " - All Smoothing Methods and Residuals'" << std::endl;
    overviewScript << std::endl;
    
    overviewScript << "set title 'All Smoothing Methods Comparison' font 'Arial,16'" << std::endl;
    overviewScript << "set xlabel 'Time Index' font 'Arial,14'" << std::endl;
    overviewScript << "set ylabel 'Price' font 'Arial,14'" << std::endl;
    overviewScript << "set grid" << std::endl;
    overviewScript << "set key outside right top vertical box" << std::endl;
    overviewScript << "set key spacing 2" << std::endl;
    overviewScript << "set key width -15" << std::endl;
    overviewScript << "set key font ',12'" << std::endl;
    overviewScript << "set xrange [*:*]" << std::endl;
    overviewScript << "set yrange [*:*]" << std::endl;
    overviewScript << std::endl;
    
    overviewScript << "set style line 1 lw 1 lc rgb '#000000' dt 1" << std::endl;  // Original - black thin
    overviewScript << "set style line 2 lw 1 lc rgb '#FF0000' dt 1" << std::endl;  // MA5 - red thin
    overviewScript << "set style line 3 lw 1 lc rgb '#00AA00' dt 1" << std::endl;  // WMA5 - green thin
    overviewScript << "set style line 4 lw 1 lc rgb '#0000FF' dt 1" << std::endl;  // WMA7 - blue thin
    overviewScript << "set style line 5 lw 1 lc rgb '#FF00FF' dt 1" << std::endl;  // CA12 - magenta thin
    overviewScript << "set style line 6 lw 1 lc rgb '#FF8C00' dt 1" << std::endl;  // Exp - orange thin
    overviewScript << std::endl;
    
    overviewScript << "plot '" << originalFile << "' using 1:3 with lines ls 1 title 'Original', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[0].getName()) << ".dat' using 1:3 with lines ls 2 title 'MA5', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[1].getName()) << ".dat' using 1:3 with lines ls 3 title 'WMA5', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[2].getName()) << ".dat' using 1:3 with lines ls 4 title 'WMA7', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[3].getName()) << ".dat' using 1:3 with lines ls 5 title 'CA12', \\" << std::endl;
    overviewScript << "     'plots/series/" << sanitizeFilename(smoothedSeries[4].getName()) << ".dat' using 1:3 with lines ls 6 title 'Exp'" << std::endl;
    overviewScript << std::endl;
    
    overviewScript << "set title 'Residuals (Original - Smoothed)' font 'Arial,16'" << std::endl;
    overviewScript << "set xlabel 'Time Index' font 'Arial,14'" << std::endl;
    overviewScript << "set ylabel 'Residual Value' font 'Arial,14'" << std::endl;
    overviewScript << "set grid" << std::endl;
    overviewScript << "set key outside right top vertical box" << std::endl;
    overviewScript << "set key spacing 2" << std::endl;
    overviewScript << "set key width -15" << std::endl;
    overviewScript << "set key font ',12'" << std::endl;
    overviewScript << "set xrange [GPVAL_X_MIN:GPVAL_X_MAX]" << std::endl; // Match x-range with first plot
    overviewScript << std::endl;
    
    overviewScript << "plot 'plots/series/residual_" << sanitizeFilename(smoothedSeries[0].getName()) << ".dat' using 1:3 with lines ls 2 lw 0.5 title 'MA5 Residuals', \\" << std::endl;
    overviewScript << "     'plots/series/residual_" << sanitizeFilename(smoothedSeries[1].getName()) << ".dat' using 1:3 with lines ls 3 lw 0.5 title 'WMA5 Residuals', \\" << std::endl;
    overviewScript << "     'plots/series/residual_" << sanitizeFilename(smoothedSeries[2].getName()) << ".dat' using 1:3 with lines ls 4 lw 0.5 title 'WMA7 Residuals', \\" << std::endl;
    overviewScript << "     'plots/series/residual_" << sanitizeFilename(smoothedSeries[3].getName()) << ".dat' using 1:3 with lines ls 5 lw 0.5 title 'CA12 Residuals', \\" << std::endl;
    overviewScript << "     'plots/series/residual_" << sanitizeFilename(smoothedSeries[4].getName()) << ".dat' using 1:3 with lines ls 6 lw 0.5 title 'Exp Residuals'" << std::endl;
    
    overviewScript << "unset multiplot" << std::endl;
    overviewScript.close();
    
    std::string overviewCommand = "gnuplot \"" + overviewScriptFile + "\"";
    int overviewResult = std::system(overviewCommand.c_str());
    
    if (overviewResult == 0) {
        std::cout << "Overview plot with residuals saved as: plots/series/smoothing_overview.png" << std::endl;
    } else {
        std::cerr << "Error: GNU Plot execution failed for overview plot" << std::endl;
    }
    
    std::string individualResidualsScriptFile = "plots/series/individual_residuals_plot.gnu";
    std::ofstream individualResidualsScript(individualResidualsScriptFile);
    individualResidualsScript << "set terminal pngcairo size 6400,3600 enhanced font 'Arial,14'" << std::endl;
    individualResidualsScript << "set output 'plots/series/individual_residuals.png'" << std::endl;
    individualResidualsScript << "set multiplot layout 3,2 title 'Individual Residuals Analysis: " << series.getName() << "'" << std::endl;
    individualResidualsScript << std::endl;
    
    for (size_t i = 0; i < smoothedSeries.size(); ++i) {
        individualResidualsScript << "set title '" << smoothedSeries[i].getName() << " Residuals' font 'Arial,16'" << std::endl;
        individualResidualsScript << "set xlabel 'Time Index' font 'Arial,14'" << std::endl;
        individualResidualsScript << "set ylabel 'Residual Value' font 'Arial,14'" << std::endl;
        individualResidualsScript << "set grid" << std::endl;
        individualResidualsScript << "plot 'plots/series/residual_" << sanitizeFilename(smoothedSeries[i].getName()) << ".dat' using 1:3 with lines lw 1 lc rgb '" 
                                  << getColorByIndex(i) << "' title 'Residuals'" << std::endl;
        individualResidualsScript << std::endl;
    }
    
    individualResidualsScript << "unset multiplot" << std::endl;
    individualResidualsScript.close();
    
    std::string individualResidualsCommand = "gnuplot \"" + individualResidualsScriptFile + "\"";
    int individualResult = std::system(individualResidualsCommand.c_str());
    
    if (individualResult == 0) {
        std::cout << "Individual residuals plot saved as: plots/series/individual_residuals.png" << std::endl;
    } else {
        std::cerr << "Error: GNU Plot execution failed for individual residuals plot" << std::endl;
    }
    
    cleanupDataFiles(dataFiles);
    cleanupDataFiles(residualFiles);
    if (std::filesystem::exists(scriptFile)) {
        std::filesystem::remove(scriptFile);
    }
    if (std::filesystem::exists(overviewScriptFile)) {
        std::filesystem::remove(overviewScriptFile);
    }
    if (std::filesystem::exists(individualResidualsScriptFile)) {
        std::filesystem::remove(individualResidualsScriptFile);
    }
}

std::string SeriesControl::getColorByIndex(size_t index) const {
    switch (index) {
        case 0: return "#FF0000";
        case 1: return "#00AA00";
        case 2: return "#0000FF"; 
        case 3: return "#FF00FF";
        case 4: return "#FF8C00"; 
        default: return "#000000";
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