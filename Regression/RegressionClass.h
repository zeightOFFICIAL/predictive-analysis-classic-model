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

struct MultipleRegressionMetrics {
    std::vector<double> coefficients;
    double R2;
    double adjustedR2;
    double Fstatistic;
    double FpValue;
    std::vector<double> tStatistics;
    std::vector<double> tpValues;
    std::vector<double> standardErrors;
    std::vector<double> residuals;
    std::vector<double> fittedValues;
    double RSS;
    double TSS;
    double ESS;
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
    static std::string sanitizeFilename(const std::string& name);

    static MultipleRegressionMetrics calculateMultipleRegression(
        const std::vector<std::vector<double>>& predictors,
        const std::vector<double>& response);
    
    static void printMultipleRegressionResults(const MultipleRegressionMetrics& results,
                                             const std::vector<std::string>& predictorNames);
    
    static double calculatePValue(double statistic, int df);
    static double calculateFStatistic(double ESS, double RSS, int p, int n);
    static std::string getSignificanceStars(double pValue);
    static void calculateStandardErrors(MultipleRegressionMetrics& results,
                                            const std::vector<std::vector<double>>& X,
                                            const std::vector<double>& y);
    static double calculatePValue(double statistic, int df1, int df2);
};