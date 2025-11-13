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

void SeriesControl::analyzeTrends() const {
    std::cout << "=== TREND ANALYSIS ===" << std::endl;
    std::cout << "Series: " << series.getName() << std::endl;
    std::cout << "Number of points: " << series.size() << std::endl;
    
    std::cout << "\n1. MEAN DIFFERENCES METHOD:" << std::endl;
    auto [has_trend_mean, t_statistic] = series.checkTrendMeanDifferences();
    
    std::cout << "Trend detected: " << (has_trend_mean ? "YES" : "NO") << std::endl;
    
    if (has_trend_mean) {
        std::cout << "Trend direction: " << (t_statistic > 0 ? "DECREASING" : "INCREASING") << std::endl;
    }
    
    std::cout << "\n2. FOSTER-STEWART METHOD:" << std::endl;
    auto [has_trend_foster, foster_statistic] = series.checkTrendFosterStewart();
    
    std::cout << "Trend detected: " << (has_trend_foster ? "YES" : "NO") << std::endl;
    
    std::cout << "\n=== FINAL CONCLUSION ===" << std::endl;
    if (has_trend_mean && has_trend_foster) {
        std::cout << "STRONG EVIDENCE of trend presence" << std::endl;
    } else if (has_trend_mean || has_trend_foster) {
        std::cout << "WEAK EVIDENCE of trend presence" << std::endl;
    } else {
        std::cout << "NO EVIDENCE of trend presence" << std::endl;
    }
    std::cout << "=================================" << std::endl << std::endl;
}

void SeriesControl::plotDecomposition(int period) const {
    std::cout << "=== TIME SERIES DECOMPOSITION ===" << std::endl;
    std::cout << "Series: " << series.getName() << std::endl;
    std::cout << "Number of points: " << series.size() << std::endl;
    std::cout << "Period: " << period << std::endl;
    
    auto decomposition = series.decomposeTimeSeries(period);
    
    std::filesystem::create_directories("plots/decomposition");
    
    // Сохраняем данные для построения графиков
    std::vector<std::string> dataFiles;
    
    // 1. Данные для трендов - ВСЕ ЛИНИИ СПЛОШНЫЕ
    std::string trendFile = "plots/decomposition/trend.dat";
    std::ofstream trend(trendFile);
    trend << "# Index Original PreliminaryTrend PrimaryTrend FinalTrend" << std::endl;
    for (size_t i = 0; i < decomposition.original.size(); ++i) {
        trend << i << " " << std::fixed << std::setprecision(6) 
              << decomposition.original[i] << " "
              << decomposition.preliminaryTrend[i] << " "
              << decomposition.primaryTrend[i] << " "
              << decomposition.finalTrend[i] << std::endl;
    }
    trend.close();
    dataFiles.push_back(trendFile);
    
    // 2. Данные для сезонных компонент
    std::string seasonalFile = "plots/decomposition/seasonal.dat";
    std::ofstream seasonal(seasonalFile);
    seasonal << "# Index FinalSeasonal" << std::endl;
    for (size_t i = 0; i < decomposition.original.size(); ++i) {
        seasonal << i << " " << std::fixed << std::setprecision(6) 
                 << decomposition.finalSeasonal[i] << std::endl;
    }
    seasonal.close();
    dataFiles.push_back(seasonalFile);
    
    // 3. Данные для остатков
    std::string residualFile = "plots/decomposition/residual.dat";
    std::ofstream residual(residualFile);
    residual << "# Index FinalResidual" << std::endl;
    for (size_t i = 0; i < decomposition.original.size(); ++i) {
        residual << i << " " << std::fixed << std::setprecision(6) 
                 << decomposition.finalResidual[i] << std::endl;
    }
    residual.close();
    dataFiles.push_back(residualFile);
    
    // Создаем скрипт для GNUplot
    std::string scriptFile = "plots/decomposition/decomposition_plot.gnu";
    std::ofstream script(scriptFile);
    
    script << "set terminal pngcairo size 7200,5400 enhanced font 'Arial,20'" << std::endl;
    script << "set output 'plots/decomposition/time_series_decomposition.png'" << std::endl;
    script << "set encoding utf8" << std::endl;
    script << "set multiplot layout 3,1 title 'Time Series Decomposition: " << series.getName() << " (Y_t = U_t + V_t + ε_t)'" << std::endl;
    script << "set lmargin 15" << std::endl;
    script << "set rmargin 15" << std::endl;
    script << "set tmargin 5" << std::endl;
    script << "set bmargin 5" << std::endl;
    script << std::endl;
    
    // ГРАФИК 1: Исходный ряд и тренды - ВСЕ ЛИНИИ СПЛОШНЫЕ
    script << "set title '1. Original Series and Trend Components' font 'Arial,28'" << std::endl;
    script << "set xlabel 'Time Index' font 'Arial,24'" << std::endl;
    script << "set ylabel 'Value' font 'Arial,24'" << std::endl;
    script << "set grid linecolor '#DDDDDD' linewidth 0.5" << std::endl;
    script << "set key outside top center horizontal font 'Arial,20'" << std::endl;
    script << "set key spacing 2" << std::endl;
    script << "set yrange [*:*]" << std::endl;
    script << "set linetype 1 lw 1.0 lc rgb '#000000'  # Original - black solid" << std::endl;
    script << "set linetype 2 lw 1.0 lc rgb '#FF0000'  # Preliminary - red solid" << std::endl;
    script << "set linetype 3 lw 1.0 lc rgb '#00AA00'  # Primary - green solid" << std::endl;
    script << "set linetype 4 lw 2.0 lc rgb '#FF8C00'  # Final - orange solid" << std::endl;
    script << "plot '" << trendFile << "' using 1:2 with lines linetype 1 title 'Original (Y_t)', \\" << std::endl;
    script << "     '" << trendFile << "' using 1:3 with lines linetype 2 title 'Preliminary Trend (MA30)', \\" << std::endl;
    script << "     '" << trendFile << "' using 1:4 with lines linetype 3 title 'Primary Trend (WMA7)', \\" << std::endl;
    script << "     '" << trendFile << "' using 1:5 with lines linetype 4 title 'Final Trend (Centered MA)'" << std::endl;
    script << std::endl;
    
    // ГРАФИК 2: Сезонные компоненты
    script << "set title '2. Seasonal Component (V_t)' font 'Arial,28'" << std::endl;
    script << "set xlabel 'Time Index' font 'Arial,24'" << std::endl;
    script << "set ylabel 'Seasonal Value' font 'Arial,24'" << std::endl;
    script << "set grid linecolor '#DDDDDD' linewidth 0.5" << std::endl;
    script << "set key outside top center horizontal font 'Arial,20'" << std::endl;
    script << "set key spacing 2" << std::endl;
    script << "set yrange [*:*]" << std::endl;
    script << "set linetype 5 lw 1.5 lc rgb '#0000FF'  # Seasonal - blue solid" << std::endl;
    script << "plot '" << seasonalFile << "' using 1:2 with lines linetype 5 title 'Final Seasonal (V_t)'" << std::endl;
    script << std::endl;
    
    // ГРАФИК 3: Остатки
    script << "set title '3. Residual Component (ε_t)' font 'Arial,28'" << std::endl;
    script << "set xlabel 'Time Index' font 'Arial,24'" << std::endl;
    script << "set ylabel 'Residual Value' font 'Arial,24'" << std::endl;
    script << "set grid linecolor '#DDDDDD' linewidth 0.5" << std::endl;
    script << "set key outside top center horizontal font 'Arial,20'" << std::endl;
    script << "set key spacing 2" << std::endl;
    script << "set yrange [*:*]" << std::endl;
    script << "set linetype 6 lw 1.0 lc rgb '#8B4513'  # Residual - brown solid" << std::endl;
    script << "plot '" << residualFile << "' using 1:2 with lines linetype 6 title 'Final Residual (ε_t)'" << std::endl;
    script << std::endl;
    
    script << "unset multiplot" << std::endl;
    script.close();
    
    // Запускаем GNUplot
    std::string command = "gnuplot \"" + scriptFile + "\"";
    int result = std::system(command.c_str());
    
    if (result == 0) {
        std::cout << "Decomposition plot saved as: plots/decomposition/time_series_decomposition.png" << std::endl;
        std::cout << "Resolution: 7200x5400 pixels" << std::endl;
        std::cout << "All lines are SOLID" << std::endl;
    } else {
        std::cerr << "Error: GNU Plot execution failed for decomposition plot" << std::endl;
    }
    
    // Очистка временных файлов
    cleanupDataFiles(dataFiles);
    if (std::filesystem::exists(scriptFile)) {
        std::filesystem::remove(scriptFile);
    }
    
    // Детальный анализ сезонности
    std::cout << "\n=== SEASONAL COMPONENT ANALYSIS ===" << std::endl;
    
    // Анализируем только уникальные значения сезонного паттерна
    std::vector<double> uniqueSeasonal;
    for (int i = 0; i < period; ++i) {
        uniqueSeasonal.push_back(decomposition.finalSeasonal[i]);
    }
    
    auto [minVal, maxVal] = std::minmax_element(uniqueSeasonal.begin(), uniqueSeasonal.end());
    double mean = std::accumulate(uniqueSeasonal.begin(), uniqueSeasonal.end(), 0.0) / period;
    
    std::cout << "Seasonal pattern has " << period << " unique values" << std::endl;
    std::cout << "Range: [" << *minVal << ", " << *maxVal << "]" << std::endl;
    std::cout << "Amplitude: " << (*maxVal - *minVal) << std::endl;
    std::cout << "Mean: " << mean << std::endl;
    
    // Если амплитуда слишком маленькая, возможно неправильный период
    if (std::abs(*maxVal - *minVal) < 0.1) {
        std::cout << "WARNING: Seasonal amplitude is very small. Consider adjusting period." << std::endl;
    }
    
    // Выводим первые несколько значений сезонного паттерна
    std::cout << "First 10 seasonal values: ";
    for (int i = 0; i < std::min(10, period); ++i) {
        std::cout << uniqueSeasonal[i] << " ";
    }
    std::cout << std::endl;
}