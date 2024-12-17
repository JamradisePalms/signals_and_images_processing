#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace std;
const double m_pi = 3.141592;

vector<double> convolution_1d(const vector<double>& signal, const vector<double>& kernel) {
    int signal_length = signal.size();
    int kernel_length = kernel.size();
    int output_length = signal_length + kernel_length - 1;
    
    vector<double> output(output_length, 0.0);
    
    vector<double> padded_signal(signal_length + 2 * (kernel_length - 1), 0.0);
    
    for (int i = 0; i < signal_length; i++) {
        padded_signal[i + kernel_length - 1] = signal[i];
    }
    
    for (int i = 0; i < output_length; i++) {
        for (int j = 0; j < kernel_length; j++) {
            output[i] += padded_signal[i + j] * kernel[kernel_length - 1 - j];
        }
    }
    
    return output;
}

vector<double> sin_signal(double F, double T, double phi, int N) {
    vector<double> signal;
    double dt = 1 / F;
    for (int i = 0; i < F * N; i++) {
        double t = i * dt;
        signal.push_back(sin(2 * m_pi * t / T + F * phi));
    }
    return signal;
}

vector<double> triangle_signal(double F, double T, double phi, double t, int N) {
    vector<double> signal;
    double dt = 1 / F;
    int i = 0;
    int period = T * F + 1;
    int flag = 0;
    while (i < N * F) {
        double time = i * dt;
        if (period > T * F && flag == 0) {
            period = 0;
            flag = 1;
            for (int j = 0; j < t * F; j++ ) {
                if (j < t * F / 2) {
                    signal.push_back(j * 2.0 / (t * F));
                    i++;
                    period++;
                }
                else {
                    signal.push_back(j * (-2.0 / (t * F)) + 2.0);
                    i++;
                    period++;
                }
            }
        } else {
            signal.push_back(0.0);
            i++;
            period++;
            if (period >= T * F) {
                flag = 0;
            }
        }
    }

    return signal;
}

void save_to_csv_double(const vector<double>& signal, double F, int N, const string& filename) {
    std::ofstream csv_file;
    csv_file.open(filename);
    csv_file << "value\n";
    for (int i = 0; i < F * N; i++) {
        csv_file << signal[i] << "\n";
    }
    csv_file.close();
}

vector<double> low_pass_filter(const vector<double>& signal, int window_size) {
    vector<double> kernel(window_size, 1.0 / window_size);
    return convolution_1d(signal, kernel);
}

vector<double> high_pass_filter(const vector<double>& signal) {
    vector<double> kernel = { -1, 2, -1 };
    // vector<double> kernel = { -1, -2, 0, 2, 1 };
    return convolution_1d(signal, kernel);
}

vector<double> impulse_noise(const vector<double>& signal, double P) {
    vector<double> noisy_signal = signal;

    double max_value = *max_element(signal.begin(), signal.end());
    double min_value = *min_element(signal.begin(), signal.end());

    for (int i = 0; i < signal.size(); i++) {
        double rand_prob = static_cast<double>(rand()) / RAND_MAX;
        if(rand_prob > P) {
            double prob = static_cast<double>(rand()) / RAND_MAX;
            if (prob > 0.5) noisy_signal[i] = max_value;
            else noisy_signal[i] = min_value;
        }
    }
    return noisy_signal;
}

vector<double> median_filter_1d(const vector<double>& signal, int aperture_size) {
    int n = signal.size();
    vector<double> filtered_signal(n, 0.0);

    int half_aperture = aperture_size / 2;

    for (int i = 0; i < n; ++i) {
        int start = max(i - half_aperture, 0);
        int end = min(i + half_aperture, n - 1);

        vector<double> window(signal.begin() + start, signal.begin() + end + 1);

        sort(window.begin(), window.end());
        int window_size = window.size();
        if (window_size % 2 == 1) {
            filtered_signal[i] = window[window_size / 2];
        } else {
            filtered_signal[i] = (window[window_size / 2 - 1] + window[window_size / 2]) / 2.0;
        }
    }

    return filtered_signal;
}

int main() {
    double T = 0.5;
    double phi = 0.0;
    double F = 128.0;
    double t = 0.5;
    int N = 6;
    vector<double> triangle_signal_data = triangle_signal(F, T, phi, t, N);

    vector<double> sin_signal_data = sin_signal(F, T, phi, N);

    int low_pass_window_size = 5;
    vector<double> low_passed_triangle = low_pass_filter(triangle_signal_data, low_pass_window_size);
    vector<double> high_passed_triangle = high_pass_filter(triangle_signal_data);

    vector<double> low_passed_sine = low_pass_filter(sin_signal_data, low_pass_window_size);
    vector<double> high_passed_sine = high_pass_filter(sin_signal_data);

    save_to_csv_double(triangle_signal_data, F, N, "original_triangle_signal.csv");
    save_to_csv_double(low_passed_triangle, F, N, "low_passed_triangle_signal.csv");
    save_to_csv_double(high_passed_triangle, F, N, "high_passed_triangle_signal.csv");

    save_to_csv_double(sin_signal_data, F, N, "original_sine_signal.csv");
    save_to_csv_double(low_passed_sine, F, N, "low_passed_sine_signal.csv");
    save_to_csv_double(high_passed_sine, F, N, "high_passed_sine_signal.csv");

    int aperture_size = 11;
    vector<double> noisy_sine_signal = impulse_noise(sin_signal_data, 0.8);
    vector<double> median_filtered_signal = median_filter_1d(noisy_sine_signal, aperture_size);
    save_to_csv_double(noisy_sine_signal, F, N, "noisy_sine_signal.csv");
    save_to_csv_double(median_filtered_signal, F, N, "filtered_sine_signal.csv");

    vector<double> noisy_triangle_signal = impulse_noise(triangle_signal_data, 0.8);
    vector<double> median_filtered_triangle_signal = median_filter_1d(noisy_triangle_signal, aperture_size);
    save_to_csv_double(noisy_triangle_signal, F, N, "noisy_triangle_signal.csv");
    save_to_csv_double(median_filtered_triangle_signal, F, N, "filtered_triangle_signal.csv");

    return 0;
}