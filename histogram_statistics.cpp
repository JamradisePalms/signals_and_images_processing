#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <algorithm>

using namespace std;

struct HistogramBin {
    double value;
    int frequency;
};

double calculateExpectedValue(const vector<HistogramBin>& histogram) {
    double sum = 0;
    int totalFrequency = 0;
    
    for (const auto& bin : histogram) {
        sum += bin.value * bin.frequency;
        totalFrequency += bin.frequency;
    }
    
    return sum / totalFrequency;
}

double calculateVariance(const vector<HistogramBin>& histogram, double expectedValue) {
    double sum = 0;
    int totalFrequency = 0;
    
    for (const auto& bin : histogram) {
        sum += pow(bin.value - expectedValue, 2) * bin.frequency;
        totalFrequency += bin.frequency;
    }
    
    return sum / totalFrequency;
}

double calculateIQR(const vector<HistogramBin>& histogram) {
    vector<double> expandedData;
    
    for (const auto& bin : histogram) {
        for (int i = 0; i < bin.frequency; i++) {
            expandedData.push_back(bin.value);
        }
    }
    
    sort(expandedData.begin(), expandedData.end());
    
    int size = expandedData.size();
    int q1Index = size / 4;
    int q3Index = (3 * size) / 4;
    
    double q1 = expandedData[q1Index];
    double q3 = expandedData[q3Index];
    
    return q3 - q1;
}

int main() {
    vector<HistogramBin> histogram;
    ifstream file("histogram.csv");
    
    string line;
    getline(file, line);
    
    while (getline(file, line)) {
        stringstream ss(line);
        string value, frequency;
        
        getline(ss, value, ',');
        getline(ss, frequency, ',');
        
        HistogramBin bin;
        bin.value = stod(value);
        bin.frequency = stoi(frequency);
        histogram.push_back(bin);
    }
    
    file.close();
    
    double expectedValue = calculateExpectedValue(histogram);
    double variance = calculateVariance(histogram, expectedValue);
    double iqr = calculateIQR(histogram);
    
    cout << "Мат.ожидание: " << expectedValue << endl;
    cout << "Дисперсия: " << variance << endl;
    cout << "Межквартильное расстояние: " << iqr << endl;
    
    return 0;
}
