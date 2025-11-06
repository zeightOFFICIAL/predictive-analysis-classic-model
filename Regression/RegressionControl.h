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

private:
    const StatisticsControl& statsController;

    void displayResults(const RegressionMetrics& results, 
                       const std::string& commodityLabel) const;
};

#endif // REGRESSIONCONTROL_H