#include "StatisticsController.h"
#include <iomanip>
#include <iostream>
#include <algorithm>

// ANSI color codes
namespace {
    const std::string RED = "\033[31m";
    const std::string GREEN = "\033[32m";
    const std::string YELLOW = "\033[33m";
    const std::string BLUE = "\033[34m";
    const std::string MAGENTA = "\033[35m";
    const std::string CYAN = "\033[36m";
    const std::string RESET = "\033[0m";
    const std::string BOLD = "\033[1m";
}

StatisticsController::StatisticsController(const StockPricesStatisticsClass& stats)
    : statsRef(stats) {}

void StatisticsController::printSectionHeader(const std::string& title) const {
    std::cout << "\n" << BOLD << BLUE << "=== " << title << " ===" << RESET << "\n";
}

void StatisticsController::printStatisticRow(const std::string& label, double value, const std::string& unit) const {
    std::cout << std::left << std::setw(20) << label 
              << ": " << GREEN << std::fixed << std::setprecision(6)
              << value << RESET << " " << unit << "\n";
}

std::string StatisticsController::interpretSkewness(double value) const {
    if (value > 0.5) return RED + "Right-skewed" + RESET;
    if (value < -0.5) return RED + "Left-skewed" + RESET;
    return GREEN + "Symmetric" + RESET;
}

std::string StatisticsController::interpretKurtosis(double value) const {
    if (value > 1.0) return RED + "Heavy-tailed" + RESET;
    if (value < -1.0) return YELLOW + "Light-tailed" + RESET;
    return GREEN + "Normal" + RESET;
}

std::string StatisticsController::formatModeOutput(const std::vector<double>& modes, bool hasSingleMode) const {
    if (modes.empty()) return YELLOW + "No mode (all unique)" + RESET;
    if (hasSingleMode) return MAGENTA + std::to_string(modes[0]) + RESET;
    
    std::string result = MAGENTA + "Multi: ";
    for (size_t i = 0; i < modes.size(); ++i) {
        if (i != 0) result += ", ";
        result += std::to_string(modes[i]);
    }
    return result + RESET;
}

void StatisticsController::showFullReport() const {
    auto commodities = statsRef.getAvailableCommodities();
    if (commodities.empty()) {
        std::cout << RED << "No statistics available!\n" << RESET;
        return;
    }

    printSectionHeader("FULL STATISTICAL REPORT");
    for (const auto& commodity : commodities) {
        showCommodityAnalysis(commodity);
    }
}

void StatisticsController::showSummaryTable() const {
    auto commodities = statsRef.getAvailableCommodities();
    if (commodities.empty()) {
        std::cout << RED << "No data available!\n" << RESET;
        return;
    }

    printSectionHeader("SUMMARY STATISTICS");
    std::cout << BOLD << std::left
              << std::setw(15) << "Commodity"
              << std::setw(12) << "Mean"
              << std::setw(12) << "Std Dev"
              << std::setw(12) << "Skewness"
              << std::setw(12) << "Kurtosis"
              << RESET << "\n";
    std::cout << std::string(63, '-') << "\n";

    for (const auto& commodity : commodities) {
        const auto& stat = statsRef.getStatistics(commodity);
        std::cout << std::left << std::setw(15) << commodity.substr(0, 14)
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << stat.mean
                  << std::setw(12) << stat.standardDeviation
                  << std::setw(12) << stat.skewness
                  << std::setw(12) << stat.kurtosis
                  << "\n";
    }
}

void StatisticsController::showCommodityAnalysis(const std::string& commodityName) const {
    const auto& stat = statsRef.getStatistics(commodityName);
    
    printSectionHeader(commodityName + " ANALYSIS");
    
    // Basic Statistics
    printStatisticRow("Data Points", stat.count);
    printStatisticRow("Mean", stat.mean);
    printStatisticRow("Median", stat.median);
    printStatisticRow("Mode", 0, formatModeOutput(stat.modes, stat.hasSingleMode));
    
    // Dispersion
    printStatisticRow("Variance", stat.variance);
    printStatisticRow("Std Deviation", stat.standardDeviation);
    printStatisticRow("Range", stat.range, "(" + std::to_string(stat.min) + " to " + std::to_string(stat.max) + ")");
    printStatisticRow("IQR", stat.iqr);
    
    // Distribution
    std::cout << std::left << std::setw(20) << "Distribution Shape" << ": "
              << "Skewness=" << stat.skewness << " (" << interpretSkewness(stat.skewness) << "), "
              << "Kurtosis=" << stat.kurtosis << " (" << interpretKurtosis(stat.kurtosis) << ")\n";
}

void StatisticsController::showComparativeAnalysis() const {
    auto commodities = statsRef.getAvailableCommodities();
    if (commodities.size() < 2) {
        std::cout << YELLOW << "Need at least 2 commodities for comparison\n" << RESET;
        return;
    }

    printSectionHeader("COMPARATIVE ANALYSIS");
    
    auto mostVolatile = *std::max_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).standardDeviation < statsRef.getStatistics(b).standardDeviation;
        });
    
    auto leastVolatile = *std::min_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).standardDeviation < statsRef.getStatistics(b).standardDeviation;
        });
    
    std::cout << "• " << BOLD << mostVolatile << RESET << " is the most volatile (σ=" 
              << statsRef.getStatistics(mostVolatile).standardDeviation << ")\n";
    std::cout << "• " << BOLD << leastVolatile << RESET << " is the least volatile (σ=" 
              << statsRef.getStatistics(leastVolatile).standardDeviation << ")\n\n";
              
    // Compare skewness
    auto mostRightSkewed = *std::max_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).skewness < statsRef.getStatistics(b).skewness;
        });
    
    std::cout << "• " << BOLD << mostRightSkewed << RESET << " has strongest right skew ("
              << statsRef.getStatistics(mostRightSkewed).skewness << ")\n";
}