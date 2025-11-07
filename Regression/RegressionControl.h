#ifndef REGRESSIONCONTROL_H
#define REGRESSIONCONTROL_H

#include "../Statistics/StatisticsControl.h"
#include "RegressionClass.h"
#include <string>

class RegressionControl {
public:
    explicit RegressionControl(const StatisticsControl& statsController);
    
    void runSilverGoldRegression() const;
    void runCopperGoldRegression() const;
    void runCommodityRegression(const std::string& commodityName, const std::string& commodityLabel) const;
    void runMultipleRegressionAllCommodities() const;
    void runMultipleRegressionSelected(const std::vector<std::string>& selectedCommodities) const;
    
private:
    const StatisticsControl& statsControl;

    void displayResults(const RegressionMetrics& results, 
                       const std::string& commodityLabel) const;
    void displayMultipleRegressionResults(const MultipleRegressionMetrics& results,
                                        const std::vector<std::string>& predictorNames) const;
    void checkMulticollinearity(const std::vector<std::vector<double>>& predictors,
                                             const std::vector<std::string>& predictorNames) const;
    void runReducedModel(const std::vector<std::vector<double>>& allPredictors,
                        const std::vector<std::string>& allPredictorNames,
                        const std::vector<double>& goldPrices,
                        const MultipleRegressionMetrics& fullResults) const;
    
    void compareModels(const MultipleRegressionMetrics& fullModel,
                      const MultipleRegressionMetrics& reducedModel,
                      size_t fullPredictors,
                      size_t reducedPredictors) const;
};

#endif