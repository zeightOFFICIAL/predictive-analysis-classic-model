#pragma once
#include <vector>
#include <string>

struct RegressionMetrics {
    double beta0 = 0.0;  
    double beta1 = 0.0;  
    double TSS = 0.0;    
    double ESS = 0.0;    
    double RSS = 0.0;   
    double R2 = 0.0;     
    std::vector<double> fittedValues;  
    std::vector<double> residuals;    
};

class RegressionClass {
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