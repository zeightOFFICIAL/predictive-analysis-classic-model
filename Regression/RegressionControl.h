#ifndef REGRESSIONCONTROL_H
#define REGRESSIONCONTROL_H

#include "../Statistics/StatisticsControl.h"
#include "RegressionClass.h"
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

#endif // REGRESSIONCONTROL_H