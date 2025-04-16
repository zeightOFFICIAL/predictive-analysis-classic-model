#ifndef REGRESSIONCONTROLLER_H
#define REGRESSIONCONTROLLER_H

#include "StatisticsController.h"
#include "../Regression/RegressionAnalysis.h"
#include <string>

class RegressionController {
public:
    explicit RegressionController(const StatisticsController& statsController);
    
    void runSilverGoldRegression() const;
    void runCopperGoldRegression() const;
    void runCommodityRegression(const std::string& commodityName, const std::string& commodityLabel) const;

private:
    const StatisticsController& statsController;

    void displayResults(const RegressionMetrics& results, 
                       const std::string& commodityLabel) const;
};

#endif // REGRESSIONCONTROLLER_H