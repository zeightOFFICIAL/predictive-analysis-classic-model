#pragma once
#include <vector>
#include <string>

struct RegressionMetrics {
    double beta0 = 0.0;  // Intercept
    double beta1 = 0.0;  // Slope
    double TSS = 0.0;    // Total Sum of Squares
    double ESS = 0.0;    // Explained Sum of Squares
    double RSS = 0.0;    // Residual Sum of Squares
    double R2 = 0.0;     // R-squared
    std::vector<double> fittedValues;  // Predicted values
    std::vector<double> residuals;    // Residuals (actual - predicted)
};

class RegressionAnalysis {
public:
    static RegressionMetrics calculateRegression(const std::vector<double>& X, const std::vector<double>& Y);
    static void plotResiduals(const std::vector<double>& goldPrices,
                            const std::vector<double>& otherPrices,
                            const RegressionMetrics& results,
                            const std::string& commodityName);
    static void generatePlot(const std::vector<double>& goldPrices,
                           const std::vector<double>& otherPrices,
                           const RegressionMetrics& results,
                           const std::string& commodityName);
};