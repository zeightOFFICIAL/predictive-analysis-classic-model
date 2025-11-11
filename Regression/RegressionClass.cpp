#include "RegressionClass.h"
#include <numeric>
#include <cmath>
#include <fstream>
#include <cstdio> 
#include <stdexcept>
#include <iomanip>
#include <filesystem>
#include <algorithm>
#include <iostream>
#include <string>
#include <limits>

namespace LinearAlgebra {    
    std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix) {
        if (matrix.empty()) return {};
        
        size_t rows = matrix.size();
        size_t cols = matrix[0].size();
        std::vector<std::vector<double>> result(cols, std::vector<double>(rows));
        
        for (size_t i = 0; i < rows; ++i) {
            for (size_t j = 0; j < cols; ++j) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

    std::vector<std::vector<double>> multiply(const std::vector<std::vector<double>>& A, 
                                            const std::vector<std::vector<double>>& B) {
        if (A.empty() || B.empty() || A[0].size() != B.size()) return {};
        
        size_t n = A.size();
        size_t m = B[0].size();
        size_t p = B.size();
        
        std::vector<std::vector<double>> result(n, std::vector<double>(m, 0.0));
        
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < m; ++j) {
                for (size_t k = 0; k < p; ++k) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return result;
    }

    std::vector<double> multiplyMatrixVector(const std::vector<std::vector<double>>& A, 
                                        const std::vector<double>& v) {
        if (A.empty() || v.empty() || A[0].size() != v.size()) return {};
        
        std::vector<double> result(A.size(), 0.0);
        for (size_t i = 0; i < A.size(); ++i) {
            for (size_t j = 0; j < v.size(); ++j) {
                result[i] += A[i][j] * v[j];
            }
        }
        return result;
    }

    std::vector<std::vector<double>> invertMatrix(const std::vector<std::vector<double>>& matrix) {
        size_t n = matrix.size();
        if (n == 0 || matrix[0].size() != n) return {};
        
        std::vector<std::vector<double>> aug(n, std::vector<double>(2 * n, 0.0));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                aug[i][j] = matrix[i][j];
            }
            aug[i][n + i] = 1.0;
        }
        
        for (size_t i = 0; i < n; ++i) {
            size_t maxRow = i;
            double maxVal = std::abs(aug[i][i]);
            for (size_t k = i + 1; k < n; ++k) {
                if (std::abs(aug[k][i]) > maxVal) {
                    maxVal = std::abs(aug[k][i]);
                    maxRow = k;
                }
            }
            
            if (maxVal < 1e-12) {
                return {}; 
            }
            
            if (maxRow != i) {
                std::swap(aug[i], aug[maxRow]);
            }
            
            double pivot = aug[i][i];
            for (size_t j = 0; j < 2 * n; ++j) {
                aug[i][j] /= pivot;
            }
            
            for (size_t k = 0; k < n; ++k) {
                if (k != i) {
                    double factor = aug[k][i];
                    for (size_t j = 0; j < 2 * n; ++j) {
                        aug[k][j] -= factor * aug[i][j];
                    }
                }
            }
        }
        
        std::vector<std::vector<double>> inverse(n, std::vector<double>(n));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                inverse[i][j] = aug[i][n + j];
            }
        }
        
        return inverse;
    }

    std::vector<double> solveOLS(const std::vector<std::vector<double>>& X, 
                            const std::vector<double>& y, double lambda = 0.01) {
        size_t n = X.size();
        size_t p = X[0].size();
        
        if (n == 0 || p == 0 || n != y.size()) {
            return {};
        }
        
        std::vector<std::vector<double>> XTX(p, std::vector<double>(p, 0.0));
        for (size_t i = 0; i < p; ++i) {
            for (size_t j = 0; j < p; ++j) {
                for (size_t k = 0; k < n; ++k) {
                    XTX[i][j] += X[k][i] * X[k][j];
                }
                if (i == j) {
                    XTX[i][j] += lambda;
                }
            }
        }
        
        std::vector<double> XTy(p, 0.0);
        for (size_t i = 0; i < p; ++i) {
            for (size_t k = 0; k < n; ++k) {
                XTy[i] += X[k][i] * y[k];
            }
        }
        
        auto invXTX = invertMatrix(XTX);
        if (invXTX.empty()) {
            return {};
        }
        
        std::vector<double> beta(p, 0.0);
        for (size_t i = 0; i < p; ++i) {
            for (size_t j = 0; j < p; ++j) {
                beta[i] += invXTX[i][j] * XTy[j];
            }
        }
        
        return beta;
    }

} // namespace LinearAlgebra

RegressionMetrics RegressionClass::calculateRegression(const std::vector<double>& X, const std::vector<double>& Y) {
    RegressionMetrics results;
    const size_t n = X.size();
    
    const double meanX = std::accumulate(X.begin(), X.end(), 0.0) / n;
    const double meanY = std::accumulate(Y.begin(), Y.end(), 0.0) / n;
    
    double covXY = 0.0, varX = 0.0;
    for (size_t i = 0; i < n; ++i) {
        covXY += (X[i] - meanX) * (Y[i] - meanY);
        varX += std::pow(X[i] - meanX, 2);
    }
    results.beta1 = covXY / varX;
    results.beta0 = meanY - results.beta1 * meanX;    
    results.fittedValues.resize(n);
    results.TSS = 0.0;
    results.RSS = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        results.fittedValues[i] = results.beta0 + results.beta1 * X[i];
        double residual = Y[i] - results.fittedValues[i];
        results.TSS += std::pow(Y[i] - meanY, 2);
        results.RSS += std::pow(residual, 2);
    }
    
    results.ESS = results.TSS - results.RSS;
    results.R2 = 1.0 - (results.RSS / results.TSS);
    
    return results;
}

std::string RegressionClass::sanitizeFilename(const std::string& name) {
    std::string sanitized = name;
    std::replace(sanitized.begin(), sanitized.end(), ' ', '_');
    std::transform(sanitized.begin(), sanitized.end(), sanitized.begin(), [](unsigned char c){ return std::tolower(c); });
    return sanitized;
}

void RegressionClass::plotResiduals(
    const std::vector<double>& goldPrices,
    const std::vector<double>& otherPrices,
    const RegressionMetrics& results,
    const std::string& typeName) {
    
    if (goldPrices.size() != otherPrices.size() || goldPrices.size() != results.fittedValues.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    std::string sanitized = sanitizeFilename(typeName);
    std::filesystem::create_directories("plots/regression/" + sanitized);

    const std::string residualsFile = "residuals_" + sanitized + ".dat";
    std::ofstream out(residualsFile);
    if (!out) {
        throw std::runtime_error("Failed to open residuals data file.");
    }

    for (size_t i = 0; i < goldPrices.size(); ++i) {
        double residual = otherPrices[i] - results.fittedValues[i];
        out << goldPrices[i] << " " << residual << "\n";
    }
    out.close();

    const std::string scriptFile = "plot_residuals_" + sanitized + ".gp";
    std::ofstream script(scriptFile);
    if (!script) {
        throw std::runtime_error("Failed to create Gnuplot script.");
    }

    script << "set terminal pngcairo enhanced font 'Arial,12'\n"
           << "set output 'plots/regression/" << sanitized << "/residuals.png'\n"
           << "set title 'Residuals Plot (" << typeName << " vs Gold)'\n"
           << "set xlabel 'Gold Price (USD)'\n"
           << "set ylabel 'Residuals (Actual - Predicted)'\n"
           << "set grid\n"
           << "set zeroaxis lt -1\n" 
           << "plot '" << residualsFile << "' using 1:2 with points pt 7 ps 0.5 lc rgb 'red' title 'Residuals', \\\n"
           << "     0 with lines lc rgb 'black' title 'Zero line'\n";
    script.close();

    int status = std::system(("gnuplot " + scriptFile).c_str());
    if (status != 0) {
        throw std::runtime_error("Gnuplot failed to execute.");
    }

    std::filesystem::remove(residualsFile);
    std::filesystem::remove(scriptFile);
}

std::string RegressionClass::getSignificanceStars(double pValue) {
    if (pValue < 0.001) return "***";
    if (pValue < 0.01) return "**";
    if (pValue < 0.05) return "*";
    if (pValue < 0.1) return ".";
    return "not sig";
}

void RegressionClass::generatePlot(
    const std::vector<double>& goldPrices,
    const std::vector<double>& otherPrices,
    const RegressionMetrics& results,
    const std::string& typeName) {
    
    if (goldPrices.size() != otherPrices.size() || goldPrices.size() != results.fittedValues.size()) {
        throw std::invalid_argument("Input vectors must have the same size.");
    }

    std::string sanitized = sanitizeFilename(typeName);
    std::filesystem::create_directories("plots/regression/" + sanitized);

    const std::string dataFilename = "regression_" + sanitized + ".dat";
    std::ofstream dataFile(dataFilename);
    if (!dataFile.is_open()) {
        throw std::runtime_error("Failed to open data file for writing.");
    }

    for (size_t i = 0; i < goldPrices.size(); ++i) {
        dataFile << goldPrices[i] << " " << otherPrices[i] << " " 
                 << results.fittedValues[i] << "\n";
    }
    dataFile.close();

    const std::string scriptFilename = "plot_regression_" + sanitized + ".gp";
    std::ofstream script(scriptFilename);
    if (!script.is_open()) {
        throw std::runtime_error("Failed to open Gnuplot script for writing.");
    }

    script << "set terminal pngcairo enhanced font 'Arial,12'\n"
           << "set output 'plots/regression/" << sanitized << "/regression.png'\n"
           << "set title '" << typeName << " vs Gold Regression (R² = " << std::fixed << std::setprecision(3) << results.R2 << ")'\n"
           << "set xlabel 'Gold Price (USD)'\n"
           << "set ylabel '" << typeName << " Price (USD)'\n"
           << "set grid\n"
           << "plot '" << dataFilename << "' using 1:2 with points pt 7 ps 0.5 lc rgb 'blue' title 'Actual', \\\n"
           << "     '' using 1:3 with lines lw 2 lc rgb 'red' title 'Fitted'\n";
    script.close();

    int status = std::system(("gnuplot " + scriptFilename).c_str());
    if (status != 0) {
        throw std::runtime_error("Gnuplot execution failed.");
    }

    std::filesystem::remove(dataFilename);
    std::filesystem::remove(scriptFilename);
}

void RegressionClass::calculateStandardErrors(MultipleRegressionMetrics& results,
                                            const std::vector<std::vector<double>>& X,
                                            const std::vector<double>& y) {
    size_t n = X.size();
    size_t p = X[0].size();
    
    if (n <= p || results.RSS < 1e-10) {
        return;
    }
    
    double mse = results.RSS / (n - p);
    auto XT = LinearAlgebra::transpose(X);
    auto XTX = LinearAlgebra::multiply(XT, X);
    auto invXTX = LinearAlgebra::invertMatrix(XTX);
    
    if (invXTX.empty()) {
        results.standardErrors.resize(p, std::sqrt(mse));
        results.tStatistics.resize(p, 0.0);
        results.tpValues.resize(p, 1.0);
        return;
    }
    
    results.standardErrors.resize(p);
    results.tStatistics.resize(p);
    results.tpValues.resize(p);
    
    for (size_t i = 0; i < p; ++i) {
        results.standardErrors[i] = std::sqrt(mse * std::max(invXTX[i][i], 0.0));
        
        if (results.standardErrors[i] > 1e-10) {
            results.tStatistics[i] = results.coefficients[i] / results.standardErrors[i];
            results.tpValues[i] = calculatePValue(std::abs(results.tStatistics[i]), 1, n - p);
        } else {
            results.tStatistics[i] = 0.0;
            results.tpValues[i] = 1.0;
        }
    }
}

static double incompleteBeta(double a, double b, double x) {
    const double EPS = 1e-12;
    const int MAX_ITER = 200;

    if (x <= 0.0) return 0.0;
    if (x >= 1.0) return 1.0;

    double lnBeta = std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b);

    bool use_symmetry = (x > (a + 1.0) / (a + b + 2.0));
    if (use_symmetry) x = 1.0 - x, std::swap(a, b);

    double qab = a + b;
    double qap = a + 1.0;
    double qam = a - 1.0;
    double c = 1.0;
    double d = 1.0 - qab * x / qap;
    if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
    d = 1.0 / d;
    double h = d;

    for (int m = 1; m <= MAX_ITER; ++m) {
        int m2 = 2 * m;
        double aa = m * (b - m) * x / ((qam + m2) * (a + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        c = 1.0 + aa / c;
        if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
        d = 1.0 / d;
        h *= d * c;
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
        d = 1.0 + aa * d;
        if (std::fabs(d) < std::numeric_limits<double>::min()) d = std::numeric_limits<double>::min();
        c = 1.0 + aa / c;
        if (std::fabs(c) < std::numeric_limits<double>::min()) c = std::numeric_limits<double>::min();
        d = 1.0 / d;
        double del = d * c;
        h *= del;
        if (std::fabs(del - 1.0) < EPS) break;
    }

    double front = std::exp(a * std::log(x) + b * std::log(1.0 - x) - lnBeta) / a;
    double result = front * h;
    return use_symmetry ? 1.0 - result : result;
}

double RegressionClass::calculatePValue(double statistic, int df1, int df2) {
    if (df1 <= 0 || df2 <= 0) return 1.0;
    if (statistic < 0) statistic = 0;

    double x = (df1 * statistic) / (df1 * statistic + df2);
    double a = df1 / 2.0;
    double b = df2 / 2.0;
    double pLower = incompleteBeta(a, b, x);
    double pUpper = 1.0 - pLower; 
    return pUpper;
}

void RegressionClass::printMultipleRegressionResults(
    const MultipleRegressionMetrics& results,
    const std::vector<std::string>& predictorNames) {
    
    std::cout << "\n=== MULTIPLE REGRESSION RESULTS ===\n";
    
    std::cout << "\nRegression Equation:\n";
    std::cout << "Gold Price = " << std::fixed << std::setprecision(6) << results.coefficients[0];
    for (size_t i = 1; i < results.coefficients.size(); ++i) {
        std::string sign = (results.coefficients[i] >= 0) ? " + " : " - ";
        std::string predName = (i - 1 < predictorNames.size()) ? predictorNames[i - 1] : "X" + std::to_string(i);
        std::cout << sign << std::abs(results.coefficients[i]) << " x " << predName;
    }
    std::cout << "\n";
    
    std::cout << "\nModel Quality Metrics:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "R-squared: " << results.R2 << " (" << (results.R2 * 100) << "%)\n";
    std::cout << "Adjusted R-squared: " << results.adjustedR2 << "\n";
    std::cout << "F-statistic: " << results.Fstatistic << "\n";
    std::cout << "Total Observations: " << results.fittedValues.size() << "\n";
    std::cout << "Number of Predictors: " << (results.coefficients.size() - 1) << "\n";
    
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "Total Sum of Squares (TSS): " << results.TSS << "\n";
    std::cout << "Explained Sum of Squares (ESS): " << results.ESS << "\n";
    std::cout << "Residual Sum of Squares (RSS): " << results.RSS << "\n";
    
    std::cout << "\nCoefficient Significance:\n";
    std::cout << std::left << std::setw(15) << "Variable" 
              << std::setw(15) << "Coefficient" 
              << std::setw(12) << "Std Error" 
              << std::setw(10) << "t-stat" 
              << std::setw(12) << "p-value" 
              << "Significance\n";
    
    std::cout << std::string(70, '-') << "\n";
    
    std::string sig = getSignificanceStars(results.tpValues[0]);
    std::cout << std::left << std::setw(15) << "Intercept"
              << std::fixed << std::setprecision(6) << std::setw(15) << results.coefficients[0]
              << std::setw(12) << results.standardErrors[0]
              << std::setw(10) << results.tStatistics[0]
              << std::setw(12) << results.tpValues[0]
              << sig << "\n";
    
    for (size_t i = 1; i < results.coefficients.size(); ++i) {
        std::string varName = (i - 1 < predictorNames.size()) ? 
                            predictorNames[i - 1] : "Predictor" + std::to_string(i);
        sig = getSignificanceStars(results.tpValues[i]);
        
        std::cout << std::left << std::setw(15) << varName
                  << std::setw(15) << results.coefficients[i]
                  << std::setw(12) << results.standardErrors[i]
                  << std::setw(10) << results.tStatistics[i]
                  << std::setw(12) << results.tpValues[i]
                  << sig << "\n";
    }
}

MultipleRegressionMetrics RegressionClass::calculateMultipleRegression(
    const std::vector<std::vector<double>>& predictors,
    const std::vector<double>& response) {
    
    MultipleRegressionMetrics results;
    
    if (predictors.empty() || response.empty()) {
        return results;
    }
    
    size_t n = response.size();
    size_t p = predictors.size();
    
    for (const auto& predictor : predictors) {
        if (predictor.size() != n) {
            return results;
        }
    }
    
    std::vector<std::vector<double>> scaledPredictors = predictors;
    std::vector<double> means(p, 0.0);
    std::vector<double> stddevs(p, 1.0);
    
    for (size_t j = 0; j < p; ++j) {
        double mean = 0.0;
        for (size_t i = 0; i < n; ++i) {
            mean += scaledPredictors[j][i];
        }
        means[j] = mean / n;
        
        double variance = 0.0;
        for (size_t i = 0; i < n; ++i) {
            variance += (scaledPredictors[j][i] - means[j]) * (scaledPredictors[j][i] - means[j]);
        }
        stddevs[j] = std::sqrt(variance / (n - 1));
        
        if (stddevs[j] > 1e-10) {
            for (size_t i = 0; i < n; ++i) {
                scaledPredictors[j][i] = (scaledPredictors[j][i] - means[j]) / stddevs[j];
            }
        }
    }
    
    std::vector<std::vector<double>> X(n, std::vector<double>(p + 1, 1.0));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < p; ++j) {
            X[i][j + 1] = scaledPredictors[j][i];
        }
    }
    
    results.coefficients = LinearAlgebra::solveOLS(X, response, 0.1);
    if (results.coefficients.empty()) {
        return results;
    }
    
    for (size_t j = 1; j <= p; ++j) {
        if (stddevs[j-1] > 1e-10) {
            results.coefficients[j] /= stddevs[j-1];
            results.coefficients[0] -= results.coefficients[j] * means[j-1];
        }
    }
    
    results.fittedValues.resize(n);
    results.residuals.resize(n);
    
    double y_mean = std::accumulate(response.begin(), response.end(), 0.0) / n;
    
    results.TSS = 0.0;
    results.RSS = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        double y_pred = results.coefficients[0];
        for (size_t j = 0; j < p; ++j) {
            y_pred += results.coefficients[j + 1] * predictors[j][i];
        }
        results.fittedValues[i] = y_pred;
        
        double residual = response[i] - y_pred;
        results.residuals[i] = residual;
        
        results.TSS += (response[i] - y_mean) * (response[i] - y_mean);
        results.RSS += residual * residual;
    }
    
    results.ESS = results.TSS - results.RSS;
    results.R2 = (results.TSS > 1e-10) ? (results.ESS / results.TSS) : 0.0;
    
    if (n > p + 1) {
        results.adjustedR2 = 1.0 - (1.0 - results.R2) * (static_cast<double>(n) - 1) / (n - p - 1);
    } else {
        results.adjustedR2 = results.R2;
    }
    
    if (p > 0 && n > p + 1 && results.RSS > 1e-10) {
        results.Fstatistic = (results.ESS / p) / (results.RSS / (n - p - 1));
        results.FpValue = calculatePValue(results.Fstatistic, p, n - p - 1);
    }
    
    calculateStandardErrors(results, X, response);
    
    return results;
}

MultipleRegressionMetrics RegressionClass::calculateMultipleRegressionWithDummy(
    const std::vector<std::vector<double>>& predictors,
    const std::vector<double>& response,
    const std::vector<int>& dummyVariable,
    bool multiplicative) {
    
    MultipleRegressionMetrics results;
    
    if (predictors.empty() || response.empty() || dummyVariable.empty()) {
        return results;
    }
    
    size_t n = response.size();
    size_t p = predictors.size();
    
    // Проверка размеров
    for (const auto& predictor : predictors) {
        if (predictor.size() != n) {
            return results;
        }
    }
    if (dummyVariable.size() != n) {
        return results;
    }
    
    // Создаем матрицу предикторов с фиктивной переменной
    std::vector<std::vector<double>> extendedPredictors;
    
    if (multiplicative) {
        // Мультипликативная модель: добавляем dummy и взаимодействия со всеми предикторами
        extendedPredictors.resize(n, std::vector<double>(2 * p + 1));
        
        for (size_t i = 0; i < n; ++i) {
            // Константа
            extendedPredictors[i][0] = 1.0;
            
            // Основные предикторы
            for (size_t j = 0; j < p; ++j) {
                extendedPredictors[i][j + 1] = predictors[j][i];
            }
            
            // Фиктивная переменная
            extendedPredictors[i][p + 1] = static_cast<double>(dummyVariable[i]);
            
            // Взаимодействия с фиктивной переменной
            for (size_t j = 0; j < p; ++j) {
                extendedPredictors[i][p + 2 + j] = predictors[j][i] * dummyVariable[i];
            }
        }
    } else {
        // Аддитивная модель: просто добавляем dummy переменную
        extendedPredictors.resize(n, std::vector<double>(p + 2));
        
        for (size_t i = 0; i < n; ++i) {
            // Константа
            extendedPredictors[i][0] = 1.0;
            
            // Основные предикторы
            for (size_t j = 0; j < p; ++j) {
                extendedPredictors[i][j + 1] = predictors[j][i];
            }
            
            // Фиктивная переменная
            extendedPredictors[i][p + 1] = static_cast<double>(dummyVariable[i]);
        }
    }
    
    // Масштабирование предикторов
    std::vector<std::vector<double>> scaledPredictors = extendedPredictors;
    size_t total_p = extendedPredictors[0].size();
    std::vector<double> means(total_p, 0.0);
    std::vector<double> stddevs(total_p, 1.0);
    
    // Не масштабируем константу и фиктивную переменную
    for (size_t j = 1; j < total_p; ++j) {
        double mean = 0.0;
        for (size_t i = 0; i < n; ++i) {
            mean += scaledPredictors[i][j];
        }
        means[j] = mean / n;
        
        double variance = 0.0;
        for (size_t i = 0; i < n; ++i) {
            variance += (scaledPredictors[i][j] - means[j]) * (scaledPredictors[i][j] - means[j]);
        }
        stddevs[j] = std::sqrt(variance / (n - 1));
        
        if (stddevs[j] > 1e-10) {
            for (size_t i = 0; i < n; ++i) {
                scaledPredictors[i][j] = (scaledPredictors[i][j] - means[j]) / stddevs[j];
            }
        }
    }
    
    // Решаем OLS
    results.coefficients = LinearAlgebra::solveOLS(scaledPredictors, response, 0.1);
    if (results.coefficients.empty()) {
        return results;
    }
    
    // Обратное масштабирование коэффициентов
    for (size_t j = 1; j < total_p; ++j) {
        if (stddevs[j] > 1e-10) {
            results.coefficients[j] /= stddevs[j];
            results.coefficients[0] -= results.coefficients[j] * means[j];
        }
    }
    
    // Расчет fitted values и остатков
    results.fittedValues.resize(n);
    results.residuals.resize(n);
    
    double y_mean = std::accumulate(response.begin(), response.end(), 0.0) / n;
    
    results.TSS = 0.0;
    results.RSS = 0.0;
    
    for (size_t i = 0; i < n; ++i) {
        double y_pred = 0.0;
        for (size_t j = 0; j < total_p; ++j) {
            y_pred += results.coefficients[j] * extendedPredictors[i][j];
        }
        results.fittedValues[i] = y_pred;
        
        double residual = response[i] - y_pred;
        results.residuals[i] = residual;
        
        results.TSS += (response[i] - y_mean) * (response[i] - y_mean);
        results.RSS += residual * residual;
    }
    
    results.ESS = results.TSS - results.RSS;
    results.R2 = (results.TSS > 1e-10) ? (results.ESS / results.TSS) : 0.0;
    
    if (n > total_p) {
        results.adjustedR2 = 1.0 - (1.0 - results.R2) * (static_cast<double>(n) - 1) / (n - total_p);
    } else {
        results.adjustedR2 = results.R2;
    }
    
    if (total_p > 1 && n > total_p && results.RSS > 1e-10) {
        results.Fstatistic = (results.ESS / (total_p - 1)) / (results.RSS / (n - total_p));
        results.FpValue = calculatePValue(results.Fstatistic, total_p - 1, n - total_p);
    }
    
    calculateStandardErrors(results, extendedPredictors, response);
    
    return results;
}

HypothesisTest RegressionClass::testDummyVariableSignificance(
    const MultipleRegressionMetrics& modelWithDummy,
    const MultipleRegressionMetrics& modelWithoutDummy,
    int n, int p_full, int p_reduced) {
    
    HypothesisTest test;
    
    test.RSS_full = modelWithDummy.RSS;
    test.RSS_reduced = modelWithoutDummy.RSS;
    test.df_full = n - p_full;
    test.df_reduced = n - p_reduced;
    
    int df_numerator = p_full - p_reduced;
    int df_denominator = n - p_full;
    
    if (df_numerator > 0 && df_denominator > 0 && test.RSS_reduced > test.RSS_full) {
        test.Fstatistic = ((test.RSS_reduced - test.RSS_full) / df_numerator) / 
                         (test.RSS_full / df_denominator);
        test.pValue = calculatePValue(test.Fstatistic, df_numerator, df_denominator);
        test.significant = (test.pValue < 0.05);
    } else {
        test.Fstatistic = 0.0;
        test.pValue = 1.0;
        test.significant = false;
    }
    
    return test;
}

void RegressionClass::printDummyVariableResults(const MultipleRegressionMetrics& results,
                                              const std::vector<std::string>& predictorNames,
                                              const HypothesisTest& hypothesisTest) {
    
    std::cout << "\n" << "=== REGRESSION WITH DUMMY VARIABLE (SANCTIONS) ===" << "\n";
    
    std::cout << "\nRegression Equation:\n";
    std::cout << "Gold Price = " << std::fixed << std::setprecision(6) << results.coefficients[0];
    
    for (size_t i = 1; i < results.coefficients.size(); ++i) {
        std::string sign = (results.coefficients[i] >= 0) ? " + " : " - ";
        
        if (i <= predictorNames.size()) {
            // Основные предикторы
            std::cout << sign << std::abs(results.coefficients[i]) << " x " << predictorNames[i - 1];
        } else if (i == predictorNames.size() + 1) {
            // Фиктивная переменная
            std::cout << sign << std::abs(results.coefficients[i]) << " x Sanctions_Dummy";
        } else {
            // Взаимодействия
            size_t interaction_idx = i - predictorNames.size() - 2;
            std::cout << sign << std::abs(results.coefficients[i]) << " x " 
                      << predictorNames[interaction_idx] << "_x_Sanctions";
        }
    }
    std::cout << "\n";
    
    std::cout << "\nModel Quality Metrics:\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "R-squared: " << results.R2 << " (" << (results.R2 * 100) << "%)\n";
    std::cout << "Adjusted R-squared: " << results.adjustedR2 << "\n";
    std::cout << "F-statistic: " << results.Fstatistic << "\n";
    
    std::cout << "\n" << "=== DUMMY VARIABLE SIGNIFICANCE TEST ===" << "\n";
    std::cout << "F-statistic for dummy variable: " << std::fixed << std::setprecision(4) 
              << hypothesisTest.Fstatistic << "\n";
    std::cout << "p-value: " << hypothesisTest.pValue << "\n";
    std::cout << "Significant at 5% level: " << (hypothesisTest.significant ? "YES" : "NO") << "\n";
    
    if (hypothesisTest.significant) {
        std::cout  << "✓ The dummy variable (sanctions) significantly improves the model\n";
        std::cout << "This suggests structural change in gold price relationships after 2022\n";
    } else {
        std::cout  << "∼ The dummy variable does not significantly improve the model\n";
        std::cout << "Gold price relationships remained stable across both periods\n";
    }
    
    std::cout << "\nIndividual Coefficients:\n";
    std::cout << std::left << std::setw(25) << "Variable" 
              << std::setw(15) << "Coefficient" 
              << std::setw(12) << "Std Error" 
              << std::setw(10) << "t-stat" 
              << std::setw(12) << "p-value" 
              << "Significance\n";
    
    std::cout << std::string(74, '-') << "\n";
    
    // Intercept
    std::string sig = getSignificanceStars(results.tpValues[0]);
    std::cout << std::left << std::setw(25) << "Intercept"
              << std::fixed << std::setprecision(6) << std::setw(15) << results.coefficients[0]
              << std::setw(12) << results.standardErrors[0]
              << std::setw(10) << results.tStatistics[0]
              << std::setw(12) << results.tpValues[0]
              << sig << "\n";
    
    // Основные предикторы
    for (size_t i = 1; i <= predictorNames.size(); ++i) {
        sig = getSignificanceStars(results.tpValues[i]);
        std::cout << std::left << std::setw(25) << predictorNames[i - 1]
                  << std::setw(15) << results.coefficients[i]
                  << std::setw(12) << results.standardErrors[i]
                  << std::setw(10) << results.tStatistics[i]
                  << std::setw(12) << results.tpValues[i]
                  << sig << "\n";
    }
    
    // Фиктивная переменная
    if (predictorNames.size() + 1 < results.coefficients.size()) {
        size_t dummy_idx = predictorNames.size() + 1;
        sig = getSignificanceStars(results.tpValues[dummy_idx]);
        std::cout << std::left << std::setw(25) << "Sanctions_Dummy"
                  << std::setw(15) << results.coefficients[dummy_idx]
                  << std::setw(12) << results.standardErrors[dummy_idx]
                  << std::setw(10) << results.tStatistics[dummy_idx]
                  << std::setw(12) << results.tpValues[dummy_idx]
                  << sig << "\n";
    }
    
    // Взаимодействия (для мультипликативной модели)
    for (size_t i = predictorNames.size() + 2; i < results.coefficients.size(); ++i) {
        size_t interaction_idx = i - predictorNames.size() - 2;
        std::string interaction_name = predictorNames[interaction_idx] + "_x_Sanctions";
        sig = getSignificanceStars(results.tpValues[i]);
        std::cout << std::left << std::setw(25) << interaction_name
                  << std::setw(15) << results.coefficients[i]
                  << std::setw(12) << results.standardErrors[i]
                  << std::setw(10) << results.tStatistics[i]
                  << std::setw(12) << results.tpValues[i]
                  << sig << "\n";
    }
}