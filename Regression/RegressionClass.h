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

struct HypothesisTest {
    double Fstatistic;
    double pValue;
    double RSS_full;
    double RSS_reduced;
    int df_full;
    int df_reduced;
    bool significant;
};

struct WaldTest {
    double chi2Statistic;
    double pValue;
    bool significant;
    std::vector<int> testedIndices;
    std::string hypothesis;
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
    
    static std::string getSignificanceStars(double pValue);
    static void calculateStandardErrors(MultipleRegressionMetrics& results,
                                            const std::vector<std::vector<double>>& X,
                                            const std::vector<double>& y);
    static double calculatePValue(double statistic, int df1, int df2);

    static MultipleRegressionMetrics calculateMultipleRegressionWithDummy(
        const std::vector<std::vector<double>>& predictors,
        const std::vector<double>& response,
        const std::vector<int>& dummyVariable,
        bool multiplicative = false);
    
    static HypothesisTest testDummyVariableSignificance(
        const MultipleRegressionMetrics& modelWithDummy,
        const MultipleRegressionMetrics& modelWithoutDummy,
        int n, int p_full, int p_reduced);
    
    static void printDummyVariableResults(const MultipleRegressionMetrics& results,
                                        const std::vector<std::string>& predictorNames,
                                        const HypothesisTest& hypothesisTest);

    static WaldTest performWaldTest(
        const MultipleRegressionMetrics& model,
        const std::vector<int>& coefficientIndices,
        const std::vector<double>& hypothesizedValues = {});
    
    static WaldTest performWaldTestForDummy(
        const MultipleRegressionMetrics& modelWithDummy,
        const std::vector<std::string>& predictorNames,
        bool multiplicative = false);
    
    static void printWaldTestResults(const WaldTest& waldTest);
};