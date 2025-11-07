#ifndef RECORDCLASS_H
#define RECORDCLASS_H

#include <string>
#include <vector>
#include <map>
#include <ctime>

class RecordClass {
public:
    static const std::string WTI_OIL;
    static const std::string GOLD;
    static const std::string SILVER;
    static const std::string NATURAL_GAS;
    static const std::string CORN;
    static const std::string WHEAT;
    static const std::string SOYBEAN;
    static const std::string COPPER;
    static const std::string PLATINUM;
    static const std::string PALLADIUM;

    bool loadFromCSV(const std::string& filename);
    float getPrice(const std::string& commodity, const std::string& date) const;
    float getPrice(const std::string& commodity, time_t date) const;
    std::map<time_t, float> getAllPrices(const std::string& commodity) const;
    
    std::vector<std::string> getAllCommodities() const;
    std::vector<time_t> getAllDates() const;
    time_t getEarliestDate() const;
    time_t getLatestDate() const;
    std::string formatDate(time_t date) const;

private:
    std::map<std::string, std::map<time_t, float>> priceData;
    std::vector<time_t> dates;

    time_t parseDate(const std::string& dateStr) const;
};

#endif