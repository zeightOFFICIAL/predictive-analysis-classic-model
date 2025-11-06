#include "StatisticsControl.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <math.h>
#include <filesystem>

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

StatisticsControl::StatisticsControl(const StatisticsClass& stats)
    : statsRef(stats) {}

void StatisticsControl::printSectionHeader(const std::string& title) const {
    std::cout << "\n" << BOLD << BLUE << "=== " << title << " ===" << RESET << "\n";
}

void StatisticsControl::printStatisticRow(const std::string& label, double value, const std::string& unit) const {
    std::cout << std::left << std::setw(20) << label 
              << ": " << GREEN << std::fixed << std::setprecision(6)
              << value << RESET << " " << unit << "\n";
}

std::string StatisticsControl::interpretSkewness(double value) const {
    if (value > 0.5) return RED + "Right-skewed" + RESET;
    if (value < -0.5) return RED + "Left-skewed" + RESET;
    return GREEN + "Symmetric" + RESET;
}

std::string StatisticsControl::interpretKurtosis(double value) const {
    if (value > 1.0) return RED + "Heavy-tailed" + RESET;
    if (value < -1.0) return YELLOW + "Light-tailed" + RESET;
    return GREEN + "Normal" + RESET;
}

std::string StatisticsControl::formatModeOutput(const std::vector<double>& modes, bool hasSingleMode) const {
    if (modes.empty()) return YELLOW + "No mode (all unique)" + RESET;
    if (hasSingleMode) return MAGENTA + std::to_string(modes[0]) + RESET;
    
    std::string result = MAGENTA + "Multi: ";
    for (size_t i = 0; i < modes.size(); ++i) {
        if (i != 0) result += ", ";
        result += std::to_string(modes[i]);
    }
    return result + RESET;
}

void StatisticsControl::showFullReport() const {
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

void StatisticsControl::showSummaryTable() const {
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

void StatisticsControl::showCommodityAnalysis(const std::string& commodityName) const {
    if (!commodityExists(commodityName)) {
        std::cout << RED << "Commodity not found: " << commodityName << RESET << "\n";
        return;
    }

    const auto& stat = statsRef.getStatistics(commodityName);
    auto getCommodityName = [](const std::string& code) -> std::string {
        if (code == RecordClass::WTI_OIL) return "Crude Oil";
        if (code == RecordClass::GOLD) return "Gold";
        if (code == RecordClass::SILVER) return "Silver";
        if (code == RecordClass::NATURAL_GAS) return "Natural Gas";
        if (code == RecordClass::CORN) return "Corn";
        if (code == RecordClass::WHEAT) return "Wheat";
        if (code == RecordClass::SOYBEAN) return "Soybean";
        if (code == RecordClass::COPPER) return "Copper";
        if (code == RecordClass::PLATINUM) return "Platinum";
        if (code == RecordClass::PALLADIUM) return "Palladium";
        return "N/A";
    };
    const std::string commodity = getCommodityName(commodityName);
    printSectionHeader(commodity + " Analysis");
    
    printStatisticRow("Data Points", stat.count);
    printStatisticRow("Mean", stat.mean);
    printStatisticRow("Median", stat.median);
    printStatisticRow("Mode", 0, formatModeOutput(stat.modes, stat.hasSingleMode));
    
    printStatisticRow("Variance", stat.variance);
    printStatisticRow("Std Deviation", stat.standardDeviation);
    printStatisticRow("Range", stat.range, "(" + std::to_string(stat.min) + " to " + std::to_string(stat.max) + ")");
    printStatisticRow("IQR", stat.iqr);
    
    std::cout << std::left << std::setw(20) << "Distribution Shape" << ": "
              << "Skewness= " << YELLOW << stat.skewness << RESET << " (" << interpretSkewness(stat.skewness) << "), "
              << "Kurtosis= " << YELLOW << stat.kurtosis << RESET << " (" << interpretKurtosis(stat.kurtosis) << ")\n";

    std::cout << std::left << std::setw(20) << "Shapiro-Wilk" << ": "
              << "Shapiro W= " << YELLOW << stat.shapiroW << RESET << ", "
              << "Shapiro p= " << YELLOW << stat.shapiroP << RESET << " (" << stat.normality << ")\n";
}

void StatisticsControl::showComparativeAnalysis() const {
    auto commodities = statsRef.getAvailableCommodities();
    if (commodities.size() < 2) {
        std::cout << YELLOW << "Need at least 2 commodities for comparison\n" << RESET;
        return;
    }

    printSectionHeader("MARKET VOLATILITY AND SKEW ANALYSIS");
    
    auto getCommodityName = [](const std::string& code) {
        if (code == RecordClass::WTI_OIL) return "Crude Oil";
        if (code == RecordClass::GOLD) return "Gold";
        if (code == RecordClass::SILVER) return "Silver";
        if (code == RecordClass::NATURAL_GAS) return "Natural Gas";
        if (code == RecordClass::CORN) return "Corn";
        if (code == RecordClass::WHEAT) return "Wheat";
        if (code == RecordClass::SOYBEAN) return "Soybean";
        if (code == RecordClass::COPPER) return "Copper";
        if (code == RecordClass::PLATINUM) return "Platinum";
        if (code == RecordClass::PALLADIUM) return "Palladium";
        return "N/A";
    };

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
              
    auto mostRightSkewed = *std::max_element(commodities.begin(), commodities.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).skewness < statsRef.getStatistics(b).skewness;
        });
    
    std::cout << "- " << BOLD << getCommodityName(mostRightSkewed) << RESET << " has strongest right skew ("
              << std::fixed << std::setprecision(4)
              << statsRef.getStatistics(mostRightSkewed).skewness << ")\n";
}

std::vector<std::string> StatisticsControl::getAvailableCommodities() const {
    return statsRef.getAvailableCommodities();
}

bool StatisticsControl::commodityExists(const std::string& commodityName) const {
    auto commodities = statsRef.getAvailableCommodities();
    return std::find(commodities.begin(), commodities.end(), commodityName) != commodities.end();
}

void StatisticsControl::showGoldCorrelations() const {
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
        if (commodity == RecordClass::WTI_OIL) commodityName = "Crude Oil";
        else if (commodity == RecordClass::GOLD) commodityName = "Gold";
        else if (commodity == RecordClass::SILVER) commodityName = "Silver";
        else if (commodity == RecordClass::NATURAL_GAS) commodityName = "Natural Gas";
        else if (commodity == RecordClass::CORN) commodityName = "Corn";
        else if (commodity == RecordClass::WHEAT) commodityName = "Wheat";
        else if (commodity == RecordClass::SOYBEAN) commodityName = "Soybean";
        else if (commodity == RecordClass::COPPER) commodityName = "Copper";
        else if (commodity == RecordClass::PLATINUM) commodityName = "Platinum";
        else if (commodity == RecordClass::PALLADIUM) commodityName = "Palladium";
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

std::vector<StatisticsControl::PlotData> StatisticsControl::preparePlotData() const {
    std::vector<PlotData> plots;
    
    auto goldStats = statsRef.getStatistics(RecordClass::GOLD);
    if (goldStats.count == 0) return plots;
    
    auto commodities = statsRef.getAvailableCommodities();
    
    for (const auto& commodity : commodities) {
        if (commodity == RecordClass::GOLD) continue;
        
        PlotData plot;
        plot.commodity = commodity;
        auto commodityStats = statsRef.getStatistics(commodity);
        
        if (commodityStats.count != goldStats.count) continue;
        
        auto goldPricesMap = statsRef.getDataRef().getAllPrices(RecordClass::GOLD);
        auto otherPricesMap = statsRef.getDataRef().getAllPrices(commodity);
        
        for (const auto& [date, goldPrice] : goldPricesMap) {
            auto it = otherPricesMap.find(date);
            if (it != otherPricesMap.end()) {
                plot.points.emplace_back(goldPrice, it->second);
            }
        }
        
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
    
    std::sort(plots.begin(), plots.end(), [](const auto& a, const auto& b) {
        return std::abs(a.correlation) > std::abs(b.correlation);
    });
    
    return plots;
}

void StatisticsControl::generateScatterPlotsWithGNUplot() const {
    auto plotsData = preparePlotData();
    if (plotsData.empty()) {
        std::cout << YELLOW << "No plot data available!\n" << RESET;
        return;
    }

    std::filesystem::create_directories("plots/correlation");

    exportPlotData(plotsData);
    
    for (const auto& plot : plotsData) {
        createGNUplotScript(plot);
        system("gnuplot plot_script.gp");
        std::cout << "\nGenerated plot for Gold vs " << getCommodityName(plot.commodity) 
                  << " (correlation: " << std::fixed << std::setprecision(4)
                  << plot.correlation << ")";

        std::filesystem::remove("plot_" + plot.commodity + ".dat");
    }

    std::filesystem::remove("plot_script.gp");
}

void StatisticsControl::exportPlotData(const std::vector<PlotData>& allData) const {
    for (const auto& plot : allData) {
        std::ofstream outFile("plot_" + plot.commodity + ".dat");
        for (const auto& [gold, other] : plot.points) {
            outFile << gold << " " << other << "\n";
        }
    }
}

std::string StatisticsControl::getCommodityName(const std::string& code) const {
    if (code == RecordClass::WTI_OIL) return "Crude Oil";
    if (code == RecordClass::GOLD) return "Gold";
    if (code == RecordClass::SILVER) return "Silver";
    if (code == RecordClass::NATURAL_GAS) return "Natural Gas";
    if (code == RecordClass::CORN) return "Corn";
    if (code == RecordClass::WHEAT) return "Wheat";
    if (code == RecordClass::SOYBEAN) return "Soybean";
    if (code == RecordClass::COPPER) return "Copper";
    if (code == RecordClass::PLATINUM) return "Platinum";
    if (code == RecordClass::PALLADIUM) return "Palladium";
    return code;
}

std::string StatisticsControl::getSanitizedCommodityName(const std::string& code) const {
    std::string name = getCommodityName(code);
    std::replace(name.begin(), name.end(), ' ', '_');
    return name;
}

void StatisticsControl::createGNUplotScript(const PlotData& data) const {
    
    std::ofstream script("plot_script.gp");
    
    script << "set terminal pngcairo enhanced font 'Arial,12'\n";
    script << "set output 'plots/correlation/" << getSanitizedCommodityName(data.commodity) << ".png'\n";
    script << "set title 'Gold vs " << getCommodityName(data.commodity) 
           << " (r = " << std::fixed << std::setprecision(3) << data.correlation << ")'\n";
    script << "set xlabel 'Gold Price (USD)'\n";
    script << "set ylabel '" << getCommodityName(data.commodity) << " Price'\n";
    script << "set grid\n";
    script << "plot 'plot_" << data.commodity << ".dat' with points pt 7 ps 0.5 title ''\n";
}

void StatisticsControl::generateHistograms() const {
    auto commodities = statsRef.getAvailableCommodities();
    if (commodities.empty()) {
        std::cout << YELLOW << "No data available for histograms!\n" << RESET;
        return;
    }

    std::filesystem::create_directories("plots/histograms");

    for (const auto& commodity : commodities) {
        auto prices = getCommodityPrices(commodity);
        if (prices.empty()) {
            std::cout << YELLOW << "No price data for " << commodity << RESET << "\n";
            continue;
        }

        prices.erase(std::remove_if(prices.begin(), prices.end(), 
            [](double price) { return price <= 0 || std::isnan(price); }), prices.end());
        
        if (prices.empty()) {
            std::cout << YELLOW << "No valid price data for " << commodity << RESET << "\n";
            continue;
        }

        size_t binCount = static_cast<size_t>(1 + 3.322 * log10(prices.size()));
        binCount = std::max(static_cast<size_t>(5), std::min(binCount, static_cast<size_t>(30))); // 5-30 bins
        
        auto [minIt, maxIt] = std::minmax_element(prices.begin(), prices.end());
        double minPrice = *minIt;
        double maxPrice = *maxIt;
        
        if (maxPrice - minPrice < 1e-10) {
            minPrice -= 0.5;
            maxPrice += 0.5;
        }
        
        double range = maxPrice - minPrice;
        double binWidth = range / binCount;

        std::vector<int> bins(binCount, 0);
        for (double price : prices) {
            size_t binIndex = static_cast<size_t>((price - minPrice) / binWidth);
            if (binIndex >= binCount) binIndex = binCount - 1;
            bins[binIndex]++;
        }

        std::string sanitized = getSanitizedCommodityName(commodity);
        std::string dataFile = "hist_" + sanitized + ".dat";
        std::ofstream out(dataFile);
        if (!out) {
            std::cerr << RED << "Failed to create data file for " << commodity << RESET << "\n";
            continue;
        }

        for (size_t i = 0; i < binCount; ++i) {
            double binStart = minPrice + i * binWidth;
            double binCenter = binStart + binWidth / 2.0;
            out << binCenter << " " << bins[i] << "\n";
        }
        out.close();

        std::string scriptFile = "hist_script_" + sanitized + ".gp";
        std::ofstream script(scriptFile);
        if (!script) {
            std::cerr << RED << "Failed to create Gnuplot script for " << commodity << RESET << "\n";
            std::remove(dataFile.c_str());
            continue;
        }

        script << "set terminal pngcairo enhanced font 'Arial,12' size 800,600\n"
               << "set output 'plots/histograms/" << sanitized << ".png'\n"
               << "set title 'Price Distribution: " << getCommodityName(commodity) << "'\n"
               << "set xlabel 'Price (USD)'\n"
               << "set ylabel 'Frequency'\n"
               << "set grid xtics ytics\n"
               << "set style fill solid 0.7 border rgb 'black'\n"
               << "set boxwidth " << binWidth * 0.8 << "\n"
               << "plot '" << dataFile << "' using 1:2 with boxes lc rgb 'steelblue' title 'Frequency Distribution'\n";

        script.close();

        std::string command = "gnuplot " + scriptFile;
        int status = system(command.c_str());
        if (status != 0) {
            std::cerr << RED << "Gnuplot failed for " << commodity << RESET << "\n";
        } else {
            std::cout << GREEN << "Generated histogram for " << getCommodityName(commodity) 
                      << " at plots/histograms/" << sanitized << ".png" << RESET << "\n";
        }

        std::remove(dataFile.c_str());
        std::remove(scriptFile.c_str());
    }
}

void StatisticsControl::generateBoxplots() const {
    auto commodities = statsRef.getAvailableCommodities();
    if (commodities.empty()) {
        std::cout << YELLOW << "No data available for boxplots!\n" << RESET;
        return;
    }

    std::filesystem::create_directories("plots/boxplots");

    for (const auto& commodity : commodities) {
        std::vector<double> prices = getCommodityPrices(commodity);
        if (prices.empty()) continue;

        std::string niceName   = getCommodityName(commodity);
        std::string sanitized  = getSanitizedCommodityName(commodity);

        std::sort(prices.begin(), prices.end());

        auto q = statsRef.calculateQuartiles(prices);
        double q1 = q[0], median = q[1], q3 = q[2];
        double iqr = q3 - q1;
        double lower = q1 - 1.5 * iqr;
        double upper = q3 + 1.5 * iqr;

        double whiskerLow  = *std::lower_bound(prices.begin(), prices.end(), lower);
        double whiskerHigh = *(--std::upper_bound(prices.begin(), prices.end(), upper));

        std::vector<double> outliers;
        for (double p : prices)
            if (p < lower || p > upper)
                outliers.push_back(p);

        std::string boxFile = "box_" + sanitized + ".dat";
        {
            std::ofstream f(boxFile);
            f << "0 " << whiskerLow << " " << q1 << " " << median << " " << q3 << " " << whiskerHigh << "\n";
        }

        std::string outFile = "out_" + sanitized + ".dat";
        {
            std::ofstream f(outFile);
            for (double v : outliers) f << v << "\n";
        }

        std::string scriptFile = "script_" + sanitized + ".gp";
        std::ofstream script(scriptFile);

        script << "set terminal pngcairo enhanced font 'Arial,14' size 800,600\n"
               << "set output 'plots/boxplots/" << sanitized << ".png'\n"
               << "set title 'Boxplot — " << niceName << " Prices'\n"
               << "set ylabel 'Price (USD)'\n"
               << "set grid ytics\n"
               << "unset xtics\n"
               << "set xrange [-0.6:0.6]\n"
               << "set style fill solid 0.45 border lt -1\n"
               << "set boxwidth 0.5\n"
               << "set key top left\n\n"

               << "plot '" << boxFile << "' using 1:3:2:6:5 with candlesticks lc rgb '#834ea7ff' lw 2 notitle, \\\n"
               << "     '' using 1:4:4:4:4 with candlesticks lt 1 lw 4 notitle";
        
        if (!outliers.empty()) {
            script << ", \\\n     '" << outFile << "' using (0.0 + 0.35*(rand(0)-0.5)):1 "
                   << "with points pt 7 ps 1.1 lc rgb '#e31a1c' title 'Outliers'";
        }
        script << "\n";

        script.close();

        int ret = system(("gnuplot " + scriptFile).c_str());
        if (ret == 0) {
            std::cout << GREEN << "Generated boxplot for " << getCommodityName(commodity) 
                      << " at plots/boxplots/" << sanitized << ".png" << RESET << "\n";
        } else {
            std::cerr << RED << "Gnuplot failed for " << niceName << RESET << "\n";
        }

        std::remove(boxFile.c_str());
        std::remove(outFile.c_str());
        std::remove(scriptFile.c_str());
    }
}

std::vector<double> StatisticsControl::getCommodityPrices(const std::string& commodity) const {
    auto pricesMap = statsRef.getDataRef().getAllPrices(commodity);
    std::vector<double> prices;
    for (const auto& [_, price] : pricesMap) {
        prices.push_back(price);
    }
    return prices;
}

void analyzeSilverGoldRegression(const std::vector<double>& goldPrices, 
        const std::vector<double>& silverPrices) {
    double meanGold = StatisticsClass::calculateMean(goldPrices);
    double meanSilver = StatisticsClass::calculateMean(silverPrices);

    double cov = 0.0, varGold = 0.0;
    for (size_t i = 0; i < goldPrices.size(); ++i) {
    cov += (goldPrices[i] - meanGold) * (silverPrices[i] - meanSilver);
    varGold += pow(goldPrices[i] - meanGold, 2);
    }

    double beta1 = cov / varGold; 
    double beta0 = meanSilver - beta1 * meanGold; 

    std::cout << "Actual OLS coefficients:\n";
    std::cout << "Slope (β1): " << beta1 << "\n";
    std::cout << "Intercept (β0): " << beta0 << "\n";
}