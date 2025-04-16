#include "StatisticsController.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <math.h>

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

StatisticsController::StatisticsController(const StatisticsClass& stats)
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
    if (!commodityExists(commodityName)) {
        std::cout << RED << "Commodity not found: " << commodityName << RESET << "\n";
        return;
    }

    const auto& stat = statsRef.getStatistics(commodityName);
    auto getCommodityName = [](const std::string& code) -> std::string {
        if (code == StockPricesRecordClass::WTI_OIL) return "Crude Oil";
        if (code == StockPricesRecordClass::GOLD) return "Gold";
        if (code == StockPricesRecordClass::SILVER) return "Silver";
        if (code == StockPricesRecordClass::NATURAL_GAS) return "Natural Gas";
        if (code == StockPricesRecordClass::CORN) return "Corn";
        if (code == StockPricesRecordClass::WHEAT) return "Wheat";
        if (code == StockPricesRecordClass::SOYBEAN) return "Soybean";
        if (code == StockPricesRecordClass::COPPER) return "Copper";
        if (code == StockPricesRecordClass::PLATINUM) return "Platinum";
        if (code == StockPricesRecordClass::PALLADIUM) return "Palladium";
        return "N/A";
    };
    const std::string commodity = getCommodityName(commodityName);
    printSectionHeader(commodity + "Analysis");
    
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

    printSectionHeader("MARKET VOLATILITY AND SKEW ANALYSIS");
    
    // Convert commodity codes to readable names
    auto getCommodityName = [](const std::string& code) {
        if (code == StockPricesRecordClass::WTI_OIL) return "Crude Oil";
        if (code == StockPricesRecordClass::GOLD) return "Gold";
        if (code == StockPricesRecordClass::SILVER) return "Silver";
        if (code == StockPricesRecordClass::NATURAL_GAS) return "Natural Gas";
        if (code == StockPricesRecordClass::CORN) return "Corn";
        if (code == StockPricesRecordClass::WHEAT) return "Wheat";
        if (code == StockPricesRecordClass::SOYBEAN) return "Soybean";
        if (code == StockPricesRecordClass::COPPER) return "Copper";
        if (code == StockPricesRecordClass::PLATINUM) return "Platinum";
        if (code == StockPricesRecordClass::PALLADIUM) return "Palladium";
        return "N/A";
    };

    // Compare volatility (std deviation)
    auto mostVolatile = *std::max_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).standardDeviation < statsRef.getStatistics(b).standardDeviation;
        });
    
    auto leastVolatile = *std::min_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).standardDeviation < statsRef.getStatistics(b).standardDeviation;
        });
    
    std::cout << "- " << BOLD << getCommodityName(mostVolatile) << RESET << " is the most volatile (ro = " 
              << std::fixed << std::setprecision(4) 
              << statsRef.getStatistics(mostVolatile).standardDeviation << ")\n";
    std::cout << "- " << BOLD << getCommodityName(leastVolatile) << RESET << " is the least volatile (ro = " 
              << statsRef.getStatistics(leastVolatile).standardDeviation << ")\n";
              
    // Compare skewness
    auto mostRightSkewed = *std::max_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).skewness < statsRef.getStatistics(b).skewness;
        });
    
    std::cout << "- " << BOLD << getCommodityName(mostRightSkewed) << RESET << " has strongest right skew ("
              << std::fixed << std::setprecision(4)
              << statsRef.getStatistics(mostRightSkewed).skewness << ")\n";
}

std::vector<std::string> StatisticsController::getAvailableCommodities() const {
    return statsRef.getAvailableCommodities();
}

bool StatisticsController::commodityExists(const std::string& commodityName) const {
    auto commodities = statsRef.getAvailableCommodities();
    return std::find(commodities.begin(), commodities.end(), commodityName) != commodities.end();
}

void StatisticsController::showGoldCorrelations() const {
    auto correlations = statsRef.calculateGoldCorrelations();
    
    printSectionHeader("GOLD PRICE CORRELATION ANALYSIS");
    std::cout << "Pearson correlation coefficients (r-values):\n\n";
    
    std::cout << BOLD << std::left 
              << std::setw(20) << "COMMODITY"
              << std::setw(12) << "r-VALUE"
              << "INTERPRETATION" << RESET << "\n";
    std::cout << std::string(50, '-') << "\n";
    
    for (const auto& [commodity, coeff] : correlations) {
        std::string commodityName;
        // Convert commodity codes to readable names
        if (commodity == StockPricesRecordClass::WTI_OIL) commodityName = "Crude Oil";
        else if (commodity == StockPricesRecordClass::GOLD) commodityName = "Gold";
        else if (commodity == StockPricesRecordClass::SILVER) commodityName = "Silver";
        else if (commodity == StockPricesRecordClass::NATURAL_GAS) commodityName = "Natural Gas";
        else if (commodity == StockPricesRecordClass::CORN) commodityName = "Corn";
        else if (commodity == StockPricesRecordClass::WHEAT) commodityName = "Wheat";
        else if (commodity == StockPricesRecordClass::SOYBEAN) commodityName = "Soybean";
        else if (commodity == StockPricesRecordClass::COPPER) commodityName = "Copper";
        else if (commodity == StockPricesRecordClass::PLATINUM) commodityName = "Platinum";
        else if (commodity == StockPricesRecordClass::PALLADIUM) commodityName = "Palladium";
        else commodityName = commodity;
        
        std::string interpretation;
        std::string color;
        
        if (coeff >= 0.7) {
            interpretation = "Strong positive";
            color = GREEN;
        } else if (coeff >= 0.3) {
            interpretation = "Moderate positive";
            color = CYAN;
        } else if (coeff >= 0.1) {
            interpretation = "Weak positive";
            color = YELLOW;
        } else if (coeff <= -0.7) {
            interpretation = "Strong negative";
            color = RED;
        } else if (coeff <= -0.3) {
            interpretation = "Moderate negative";
            color = MAGENTA;
        } else if (coeff <= -0.1) {
            interpretation = "Weak negative";
            color = YELLOW;
        } else {
            interpretation = "No correlation";
            color = RESET;
        }
        
        std::ostringstream valueStream;
        valueStream << std::fixed << std::setprecision(4) << coeff;
        std::string valueStr = valueStream.str();

        std::cout << std::left << std::setw(20) << commodityName
                  << std::left << std::setw(21) << color + valueStr + RESET
                  << color << interpretation << RESET << "\n";
    }
}

std::vector<StatisticsController::PlotData> StatisticsController::preparePlotData() const {
    std::vector<PlotData> plots;
    
    // Get gold prices from statistics data
    auto goldStats = statsRef.getStatistics(StockPricesRecordClass::GOLD);
    if (goldStats.count == 0) return plots;
    
    // Get all commodities except gold
    auto commodities = statsRef.getAvailableCommodities();
    
    for (const auto& commodity : commodities) {
        if (commodity == StockPricesRecordClass::GOLD) continue;
        
        PlotData plot;
        plot.commodity = commodity;
        auto commodityStats = statsRef.getStatistics(commodity);
        
        // We can't correlate if we don't have matching counts
        if (commodityStats.count != goldStats.count) continue;
        
        // Get raw prices from RecordClass (assuming we can access it)
        auto goldPricesMap = statsRef.getDataRef().getAllPrices(StockPricesRecordClass::GOLD);
        auto otherPricesMap = statsRef.getDataRef().getAllPrices(commodity);
        
        // Align prices by date
        for (const auto& [date, goldPrice] : goldPricesMap) {
            auto it = otherPricesMap.find(date);
            if (it != otherPricesMap.end()) {
                plot.points.emplace_back(goldPrice, it->second);
            }
        }
        
        // Calculate correlation if we have data
        if (!plot.points.empty()) {
            std::vector<double> x, y;
            for (const auto& [xVal, yVal] : plot.points) {
                x.push_back(xVal);
                y.push_back(yVal);
            }
            plot.correlation = StatisticsClass::calculatePearsonCorrelation(x, y);
            plots.push_back(plot);
        }
    }
    
    // Sort by absolute correlation strength
    std::sort(plots.begin(), plots.end(), [](const auto& a, const auto& b) {
        return std::abs(a.correlation) > std::abs(b.correlation);
    });
    
    return plots;
}

void StatisticsController::generateScatterPlotsWithGNUplot() const {
    auto plotsData = preparePlotData();
    exportPlotData(plotsData);
    
    for (const auto& plot : plotsData) {
        createGNUplotScript(plot);
        system("gnuplot plot_script.gp");
        std::cout << "\nGenerated plot for Gold vs " << plot.commodity 
                  << " (correlation: " << std::fixed << std::setprecision(4)
                  << plot.correlation << ")\n";
    }
}

void StatisticsController::exportPlotData(const std::vector<PlotData>& allData) const {
    for (const auto& plot : allData) {
        std::ofstream outFile("plot_" + plot.commodity + ".dat");
        for (const auto& [gold, other] : plot.points) {
            outFile << gold << " " << other << "\n";
        }
    }
}

std::string StatisticsController::getCommodityName(const std::string& code) const {
    if (code == StockPricesRecordClass::WTI_OIL) return "Crude Oil";
    if (code == StockPricesRecordClass::GOLD) return "Gold";
    if (code == StockPricesRecordClass::SILVER) return "Silver";
    if (code == StockPricesRecordClass::NATURAL_GAS) return "Natural Gas";
    if (code == StockPricesRecordClass::CORN) return "Corn";
    if (code == StockPricesRecordClass::WHEAT) return "Wheat";
    if (code == StockPricesRecordClass::SOYBEAN) return "Soybean";
    if (code == StockPricesRecordClass::COPPER) return "Copper";
    if (code == StockPricesRecordClass::PLATINUM) return "Platinum";
    if (code == StockPricesRecordClass::PALLADIUM) return "Palladium";
    return code; // Return the code if no match found
}

void StatisticsController::createGNUplotScript(const PlotData& data) const {

    
    std::ofstream script("plot_script.gp");
    
    script << "set terminal pngcairo enhanced font 'Arial,12'\n";
    script << "set output 'gold_vs_" << data.commodity << ".png'\n";
    script << "set title 'Gold vs " << getCommodityName(data.commodity) 
           << " (r = " << std::fixed << std::setprecision(3) << data.correlation << ")'\n";
    script << "set xlabel 'Gold Price (USD)'\n";
    script << "set ylabel '" << getCommodityName(data.commodity) << " Price'\n";
    script << "set grid\n";
    script << "plot 'plot_" << data.commodity << ".dat' with points pt 7 ps 0.5 title ''\n";
}

std::vector<double> StatisticsController::getCommodityPrices(const std::string& commodity) const {
    auto pricesMap = statsRef.getDataRef().getAllPrices(commodity);
    std::vector<double> prices;
    for (const auto& [_, price] : pricesMap) {
        prices.push_back(price);
    }
    return prices;
}

// After calculating correlations in StatisticsController.cpp
void analyzeSilverGoldRegression(const std::vector<double>& goldPrices, 
        const std::vector<double>& silverPrices) {
    // Calculate means
    double meanGold = StatisticsClass::calculateMean(goldPrices);
    double meanSilver = StatisticsClass::calculateMean(silverPrices);

    // Calculate covariance and variance
    double cov = 0.0, varGold = 0.0;
    for (size_t i = 0; i < goldPrices.size(); ++i) {
    cov += (goldPrices[i] - meanGold) * (silverPrices[i] - meanSilver);
    varGold += pow(goldPrices[i] - meanGold, 2);
    }

    // Calculate coefficients
    double beta1 = cov / varGold;  // Slope
    double beta0 = meanSilver - beta1 * meanGold;  // Intercept

    std::cout << "Actual OLS coefficients:\n";
    std::cout << "Slope (β1): " << beta1 << "\n";
    std::cout << "Intercept (β0): " << beta0 << "\n";
}