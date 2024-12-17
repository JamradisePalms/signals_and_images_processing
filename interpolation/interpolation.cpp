#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

vector<double> sin_signal(double F, double T, double phi, int N) {
    vector<double> signal;
    double dt = 1 / F;
    for (int i = 0; i < F * N; i++) {
        double t = i * dt;
        signal.push_back(sin(2 * 3.1415 * t / T + F * phi));
    }
    return signal;
}

vector<double> zero_order_interpolation(const vector<double>& signal, int factor) {
    vector<double> interpolated_signal;
    for (double value : signal) {
        for (int i = 0; i < factor; ++i) {
            interpolated_signal.push_back(value);
        }
    }
    return interpolated_signal;
}

vector<double> first_order_interpolation(const vector<double>& signal, int factor) {
    vector<double> interpolated_signal;
    for (size_t i = 0; i < signal.size() - 1; ++i) {
        double start = signal[i];
        double end = signal[i + 1];
        for (int j = 0; j < factor; ++j) {
            double interpolated_value = start + j * (end - start) / factor;
            interpolated_signal.push_back(interpolated_value);
        }
    }
    interpolated_signal.push_back(signal.back());
    return interpolated_signal;
}

void save_signal_to_csv(const vector<double>& signal, const string& filename) {
    ofstream file(filename);
    file << "Index,Value\n";
    for (size_t i = 0; i < signal.size(); ++i) {
        file << i << "," << signal[i] << "\n";
    }
    file.close();
}

int main() {
    double T = 0.5;
    double phi = 0.0;
    double F = 24.0;
    int N = 6;
    vector<double> sin_signal_data = sin_signal(F, T, phi, N);

    vector<double> signal = {1, 2, 3, 5, 8, 3, 5, 7, 1, 2, 10, 15, 20, 20, 9, 1};
    int interpolation_factor = 4;

    vector<double> zero_order_signal = zero_order_interpolation(signal, interpolation_factor);
    vector<double> first_order_signal = first_order_interpolation(signal, interpolation_factor);


    save_signal_to_csv(zero_order_signal, "zero_order_signal.csv");
    save_signal_to_csv(first_order_signal, "first_order_signal.csv");

    vector<double> zero_order_sine_signal = zero_order_interpolation(sin_signal_data, interpolation_factor);
    vector<double> first_order_sine_signal = first_order_interpolation(sin_signal_data, interpolation_factor);

    save_signal_to_csv(sin_signal_data, "sin_signal.csv");
    save_signal_to_csv(zero_order_sine_signal, "zero_order_sine_signal.csv");
    save_signal_to_csv(first_order_sine_signal, "first_order_sine_signal.csv");

    return 0;
}