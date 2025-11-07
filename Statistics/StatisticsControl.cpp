#include "StatisticsControl.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
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

std::string StatisticsControl::interpretShapiroWilk(double value) const {
    if (value > 0.45) return GREEN + "Normal" + RESET;
    else if (value > 0.20) return GREEN + "Borderline normal" + RESET;
    else if (value > 0.10) return YELLOW + "Moderately non-normal" + RESET;
    else return RED + "Strongly non-normal" + RESET;
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
    auto types = statsRef.getAvailableTypes();
    if (types.empty()) {
        std::cout << RED << "No statistics available!\n" << RESET;
        return;
    }

    printSectionHeader("FULL STATISTICAL REPORT");
    for (const auto& type : types) {
        showCommodityAnalysis(type);
    }
}

void StatisticsControl::showSummaryTable() const {
    auto types = statsRef.getAvailableTypes();
    if (types.empty()) {
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

    for (const auto& type : types) {
        const auto& stat = statsRef.getStatistics(type);
        std::cout << std::left << std::setw(15) << type.substr(0, 14)
                  << std::fixed << std::setprecision(4)
                  << std::setw(12) << stat.mean
                  << std::setw(12) << stat.std
                  << std::setw(12) << stat.skewness
                  << std::setw(12) << stat.kurtosis
                  << "\n";
    }
}

void StatisticsControl::showCommodityAnalysis(const std::string& typeName) const {
    if (!typeExists(typeName)) {
        std::cout << RED << "Commodity not found: " << typeName << RESET << "\n";
        return;
    }

    const auto& stat = statsRef.getStatistics(typeName);
    auto getTypeClearName = [](const std::string& code) -> std::string {
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
    
    const std::string type = getTypeClearName(typeName);
    printSectionHeader(type + " Analysis");    
    printStatisticRow("Data Points", stat.count);
    printStatisticRow("Mean", stat.mean);
    printStatisticRow("Median", stat.median);
    printStatisticRow("Mode", 0, formatModeOutput(stat.modes, stat.hasSingleMode));    
    printStatisticRow("Variance", stat.variance);
    printStatisticRow("Std Deviation", stat.std);
    printStatisticRow("Range", stat.range, "(" + std::to_string(stat.min) + " to " + std::to_string(stat.max) + ")");
    printStatisticRow("IQR", stat.iqr);    
    std::cout << std::left << std::setw(20) << "Distribution Shape" << ": "
              << "Skewness= " << YELLOW << stat.skewness << RESET << " (" << interpretSkewness(stat.skewness) << "), "
              << "Kurtosis= " << YELLOW << stat.kurtosis << RESET << " (" << interpretKurtosis(stat.kurtosis) << ")\n";
    std::cout << std::left << std::setw(20) << "Shapiro-Wilk" << ": "
              << "Shapiro W= " << YELLOW << stat.shapiroW << RESET << ", "
              << "Shapiro p= " << YELLOW << stat.shapiroP << RESET << " (" << interpretShapiroWilk(stat.shapiroP) << ")\n";
}

void StatisticsControl::showComparativeAnalysis() const {
    auto types = statsRef.getAvailableTypes();
    if (types.size() < 2) {
        std::cout << YELLOW << "Need at least 2 commodities for comparison\n" << RESET;
        return;
    }

    printSectionHeader("MARKET VOLATILITY AND SKEW ANALYSIS");    
    auto getTypeClearName = [](const std::string& code) {
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

    auto mostVolatile = *std::max_element(types.begin(), types.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).std < statsRef.getStatistics(b).std;
        });    
    auto leastVolatile = *std::min_element(types.begin(), types.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).std < statsRef.getStatistics(b).std;
        });    
    std::cout << "- " << BOLD << getTypeClearName(mostVolatile) << RESET << " is the most volatile (ro = " 
              << std::fixed << std::setprecision(4) 
              << statsRef.getStatistics(mostVolatile).std << ")\n";
    std::cout << "- " << BOLD << getTypeClearName(leastVolatile) << RESET << " is the least volatile (ro = " 
              << statsRef.getStatistics(leastVolatile).std << ")\n";              
    auto mostRightSkewed = *std::max_element(types.begin(), types.end(),
        [this](const auto& a, const auto& b) {
            return statsRef.getStatistics(a).skewness < statsRef.getStatistics(b).skewness;
        });    
    std::cout << "- " << BOLD << getTypeClearName(mostRightSkewed) << RESET << " has strongest right skew ("
              << std::fixed << std::setprecision(4)
              << statsRef.getStatistics(mostRightSkewed).skewness << ")\n";
}

std::vector<std::string> StatisticsControl::getAvailableTypes() const {
    return statsRef.getAvailableTypes();
}

bool StatisticsControl::typeExists(const std::string& typeName) const {
    auto types = statsRef.getAvailableTypes();
    return std::find(types.begin(), types.end(), typeName) != types.end();
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
    
    for (const auto& [type, coeff] : correlations) {
        std::string typeName;
        if (type == RecordClass::WTI_OIL) typeName = "Crude Oil";
        else if (type == RecordClass::GOLD) typeName = "Gold";
        else if (type == RecordClass::SILVER) typeName = "Silver";
        else if (type == RecordClass::NATURAL_GAS) typeName = "Natural Gas";
        else if (type == RecordClass::CORN) typeName = "Corn";
        else if (type == RecordClass::WHEAT) typeName = "Wheat";
        else if (type == RecordClass::SOYBEAN) typeName = "Soybean";
        else if (type == RecordClass::COPPER) typeName = "Copper";
        else if (type == RecordClass::PLATINUM) typeName = "Platinum";
        else if (type == RecordClass::PALLADIUM) typeName = "Palladium";
        else typeName = type;
        
        std::string interpretation;
        std::string color;
        
        if (coeff >= 0.7) {
            interpretation = "Strong positive";
            color = GREEN;
        } else if (coeff >= 0.4) {
            interpretation = "Moderate positive";
            color = CYAN;
        } else if (coeff >= 0.2) {
            interpretation = "Weak positive";
            color = YELLOW;
        } else if (coeff <= -0.7) {
            interpretation = "Strong negative";
            color = RED;
        } else if (coeff <= -0.4) {
            interpretation = "Moderate negative";
            color = MAGENTA;
        } else if (coeff <= -0.2) {
            interpretation = "Weak negative";
            color = YELLOW;
        } else {
            interpretation = "No correlation";
            color = RESET;
        }
        
        std::ostringstream valueStream;
        valueStream << std::fixed << std::setprecision(4) << coeff;
        std::string valueStr = valueStream.str();
        std::cout << std::left << std::setw(20) << typeName
                  << std::left << std::setw(21) << color + valueStr + RESET
                  << color << interpretation << RESET << "\n";
    }
}

std::vector<StatisticsControl::PlotData> StatisticsControl::preparePlotData() const {
    std::vector<PlotData> plots;
    
    auto goldStats = statsRef.getStatistics(RecordClass::GOLD);
    if (goldStats.count == 0) return plots;
    
    auto types = statsRef.getAvailableTypes();
    
    for (const auto& type : types) {
        if (type == RecordClass::GOLD) continue;
        
        PlotData plot;
        plot.type = type;
        auto typeStats = statsRef.getStatistics(type);
        
        if (typeStats.count != goldStats.count) continue;
        
        auto goldPricesMap = statsRef.getDataRef().getAllPrices(RecordClass::GOLD);
        auto otherPricesMap = statsRef.getDataRef().getAllPrices(type);
        
        for (const auto& [date, goldPrice] : goldPricesMap) {
            auto it = otherPricesMap.find(date);
            if (it != otherPricesMap.end()) {
                plot.points.emplace_back(static_cast<double>(goldPrice), static_cast<double>(it->second));
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
        system(("gnuplot plot_script_" + getSanitizedTypeName(plot.type) + ".gp").c_str());
        std::cout << "\nGenerated plot for Gold vs " << getTypeName(plot.type) 
                  << " (correlation: " << std::fixed << std::setprecision(4)
                  << plot.correlation << ")";

        std::filesystem::remove("plot_" + plot.type + ".dat");
        std::filesystem::remove("plot_script_" + getSanitizedTypeName(plot.type) + ".gp");
    }
}

void StatisticsControl::exportPlotData(const std::vector<PlotData>& allData) const {
    for (const auto& plot : allData) {
        std::ofstream outFile("plot_" + plot.type + ".dat");
        for (const auto& [gold, other] : plot.points) {
            outFile << gold << " " << other << "\n";
        }
    }
}

std::string StatisticsControl::getTypeName(const std::string& code) const {
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

std::string StatisticsControl::getSanitizedTypeName(const std::string& code) const {
    std::string name = getTypeName(code);
    std::replace(name.begin(), name.end(), ' ', '_');
    return name;
}

void StatisticsControl::createGNUplotScript(const PlotData& data) const {
    std::string sanitized = getSanitizedTypeName(data.type);
    std::ofstream script("plot_script_" + sanitized + ".gp");    
    script << "set terminal pngcairo enhanced font 'Arial,12'\n";
    script << "set output 'plots/correlation/" << sanitized << ".png'\n";
    script << "set title 'Gold vs " << getTypeName(data.type) 
           << " (r = " << std::fixed << std::setprecision(3) << data.correlation << ")'\n";
    script << "set xlabel 'Gold Price (USD)'\n";
    script << "set ylabel '" << getTypeName(data.type) << " Price'\n";
    script << "set grid\n";
    script << "plot 'plot_" << data.type << ".dat' with points pt 7 ps 0.5 title ''\n";
}

void StatisticsControl::generateHistograms() const {
    auto types = statsRef.getAvailableTypes();
    if (types.empty()) {
        std::cout << YELLOW << "No data available for histograms!\n" << RESET;
        return;
    }

    std::filesystem::create_directories("plots/histograms");

    for (const auto& type : types) {
        auto prices = getTypePrices(type);
        if (prices.empty()) {
            std::cout << YELLOW << "No price data for " << type << RESET << "\n";
            continue;
        }

        prices.erase(std::remove_if(prices.begin(), prices.end(), 
            [](double price) { return price <= 0 || std::isnan(price); }), prices.end());
        
        if (prices.empty()) {
            std::cout << YELLOW << "No valid price data for " << type << RESET << "\n";
            continue;
        }

        size_t binCount = static_cast<size_t>(1 + 3.322 * log10(prices.size()));
        binCount = std::max(static_cast<size_t>(5), std::min(binCount, static_cast<size_t>(30)));
        
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

        std::string sanitized = getSanitizedTypeName(type);
        std::string dataFile = "hist_" + sanitized + ".dat";
        std::ofstream out(dataFile);
        if (!out) {
            std::cerr << RED << "Failed to create data file for " << type << RESET << "\n";
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
            std::cerr << RED << "Failed to create Gnuplot script for " << type << RESET << "\n";
            std::filesystem::remove(dataFile);
            continue;
        }

        script << "set terminal pngcairo enhanced font 'Arial,12' size 800,600\n"
               << "set output 'plots/histograms/" << sanitized << ".png'\n"
               << "set title 'Price Distribution: " << getTypeName(type) << "'\n"
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
            std::cerr << RED << "Gnuplot failed for " << type << RESET << "\n";
        } else {
            std::cout << GREEN << "Generated histogram for " << getTypeName(type) 
                      << " at plots/histograms/" << sanitized << ".png" << RESET << "\n";
        }

        std::filesystem::remove(dataFile);
        std::filesystem::remove(scriptFile);
    }
}

void StatisticsControl::generateBoxplots() const {
    auto types = statsRef.getAvailableTypes();
    if (types.empty()) {
        std::cout << YELLOW << "No data available for boxplots!\n" << RESET;
        return;
    }

    std::filesystem::create_directories("plots/boxplots");

    for (const auto& type : types) {
        std::vector<double> prices = getTypePrices(type);
        if (prices.empty()) continue;

        std::string niceName   = getTypeName(type);
        std::string sanitized  = getSanitizedTypeName(type);

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
            std::cout << GREEN << "Generated boxplot for " << getTypeName(type) 
                      << " at plots/boxplots/" << sanitized << ".png" << RESET << "\n";
        } else {
            std::cerr << RED << "Gnuplot failed for " << niceName << RESET << "\n";
        }

        std::filesystem::remove(boxFile);
        std::filesystem::remove(outFile);
        std::filesystem::remove(scriptFile);
    }
}

void StatisticsControl::generateCorrelationMatrix() const {
    auto types = statsRef.getAvailableTypes();
    if (types.empty()) {
        std::cout << YELLOW << "No data available for correlation matrix!\n" << RESET;
        return;
    }

    printSectionHeader("CORRELATION MATRIX - ALL COMMODITIES");

    std::vector<std::vector<double>> allPrices;
    std::vector<std::string> typesNames;
    
    for (const auto& type : types) {
        auto prices = getTypePrices(type);
        if (!prices.empty()) {
            allPrices.push_back(prices);
            typesNames.push_back(getTypeName(type));
        }
    }

    size_t n = allPrices.size();
    if (n < 2) {
        std::cout << YELLOW << "Need at least 2 commodities with data for correlation matrix\n" << RESET;
        return;
    }

    std::vector<std::vector<double>> correlationMatrix(n, std::vector<double>(n, 0.0));

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            if (i == j) {
                correlationMatrix[i][j] = 1.0;
            } else {
                size_t minSize = std::min(allPrices[i].size(), allPrices[j].size());
                if (minSize < 2) {
                    correlationMatrix[i][j] = 0.0;
                    continue;
                }

                std::vector<double> x(minSize), y(minSize);
                std::copy(allPrices[i].begin(), allPrices[i].begin() + minSize, x.begin());
                std::copy(allPrices[j].begin(), allPrices[j].begin() + minSize, y.begin());
                
                correlationMatrix[i][j] = StatisticsClass::calculatePearsonCorrelation(x, y);
            }
        }
    }

    const int nameWidth = 15;
    const int valueWidth = 10;

    std::cout << std::left << std::setw(nameWidth) << "";
    for (size_t i = 0; i < n; ++i) {
        std::string shortName = typesNames[i].substr(0, nameWidth - 1);
        std::cout << std::setw(valueWidth) << shortName;
    }
    std::cout << "\n" << std::string(nameWidth + n * valueWidth, '-') << "\n";

    for (size_t i = 0; i < n; ++i) {
        std::cout << std::left << std::setw(nameWidth) 
                  << typesNames[i].substr(0, nameWidth - 1);
        
        for (size_t j = 0; j < n; ++j) {
            double corr = correlationMatrix[i][j];
            std::string color = RESET;
            
            if (i != j) {
                if (corr >= 0.7) color = GREEN;
                else if (corr >= 0.3) color = CYAN;
                else if (corr >= 0.1) color = YELLOW;
                else if (corr <= -0.7) color = RED;
                else if (corr <= -0.3) color = MAGENTA;
                else if (corr <= -0.1) color = YELLOW;
            }

            std::ostringstream oss;
            oss << std::fixed << std::setprecision(3) << corr;
            std::cout << color << std::setw(valueWidth) << oss.str() << RESET;
        }
        std::cout << "\n";
    }

    generateCorrelationMatrixPlot(correlationMatrix, typesNames);
}

void StatisticsControl::generateCorrelationMatrixPlot(
    const std::vector<std::vector<double>>& correlationMatrix,
    const std::vector<std::string>& typesNames) const {
    
    std::filesystem::create_directories("plots/correlation");

    std::ofstream dataFile("correlation_matrix.dat");
    if (!dataFile) {
        std::cerr << RED << "Failed to create correlation data file!\n" << RESET;
        return;
    }

    size_t n = correlationMatrix.size();
    
    dataFile << "# Correlation Matrix Data\n";
    dataFile << "# ";
    for (size_t i = 0; i < n; ++i) {
        dataFile << typesNames[i];
        if (i < n - 1) dataFile << " | ";
    }
    dataFile << "\n";

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            dataFile << std::fixed << std::setprecision(6) 
                     << correlationMatrix[i][j];
            if (j < n - 1) dataFile << " ";
        }
        dataFile << "\n";
    }
    dataFile.close();

    std::ofstream script("correlation_heatmap.gp");
    if (!script) {
        std::cerr << RED << "Failed to create GNUplot script!\n" << RESET;
        return;
    }

    script << "set terminal pngcairo enhanced font 'Arial,12' size 1000,800\n";
    script << "set output 'plots/correlation/Matrix.png'\n";
    script << "set title 'Commodity Correlation Matrix Heatmap'\n";
    script << "unset key\n";
    script << "set tic scale 0\n";
    
    script << "set palette defined (-1 'yellow', 0 'white', 1 'green')\n";
    script << "set cbrange [-1:1]\n";
    script << "set cblabel 'Correlation Coefficient'\n";
    
    script << "set view map\n";
    script << "set size ratio 1\n";
    
    script << "set xtics border in scale 0,0.001 autojustify\n";
    script << "set ytics border in scale 0,0.001 autojustify\n";
    script << "set xrange [-0.5:" << (n - 0.5) << "]\n";
    script << "set yrange [-0.5:" << (n - 0.5) << "]\n";
    
    script << "set xtics rotate by 45 right\n";
    
    script << "set xtics (";
    for (size_t i = 0; i < n; ++i) {
        script << "\"" << typesNames[i] << "\" " << i;
        if (i < n - 1) script << ", ";
    }
    script << ")\n";
    
    script << "set ytics (";
    for (size_t i = 0; i < n; ++i) {
        script << "\"" << typesNames[i] << "\" " << i;
        if (i < n - 1) script << ", ";
    }
    script << ")\n";
    
    script << "plot 'correlation_matrix.dat' matrix with image, \\\n";
    script << "     '' matrix using 1:2:($3 >= 0 ? sprintf('%.3f', $3) : '') with labels font 'Arial,8' center, \\\n";
    script << "     '' matrix using 1:2:($3 < 0 ? sprintf('%.3f', $3) : '') with labels font 'Arial,8' center\n";

    script.close();

    int result = system("gnuplot correlation_heatmap.gp");
    if (result == 0) {
        std::cout << GREEN << "Correlation matrix heatmap saved to: plots/correlation/Matrix.png" << RESET << "\n";
    } else {
        std::cerr << YELLOW << "GNUplot execution failed. Make sure GNUplot is installed.\n" << RESET;
    }

    std::filesystem::remove("correlation_matrix.dat");
    std::filesystem::remove("correlation_heatmap.gp");
}

std::vector<double> StatisticsControl::getTypePrices(const std::string& type) const {
    auto pricesMap = statsRef.getDataRef().getAllPrices(type);
    std::vector<double> prices;
    for (const auto& [_, price] : pricesMap) {
        prices.push_back(static_cast<double>(price));
    }
    return prices;
}