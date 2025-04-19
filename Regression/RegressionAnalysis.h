#ifndef REGRESSIONANALYSIS_H
#define REGRESSIONANALYSIS_H

#include <vector>
#include <string>

struct RegressionMetrics {
    double beta0;
    double beta1;
    double TSS;  // Total Sum of Squares
    double ESS;  // Explained Sum of Squares
    double RSS;  // Residual Sum of Squares
    double R2;   // R-squared
    std::vector<double> fittedValues;
};

class RegressionAnalysis {
public:
    static RegressionMetrics calculateRegression(const std::vector<double>& X, const std::vector<double>& Y);
    static void generatePlot(const std::vector<double>& goldPrices, 
                           const std::vector<double>& otherPrices,
                           const RegressionMetrics& results,
                           const std::string& commodityName);
    static void plotResiduals(const std::vector<double>& goldPrices,
                            const std::vector<double>& otherPrices,
                            const RegressionMetrics& results,
                            const std::string& commodityName);
};

#endif // REGRESSIONANALYSIS_H