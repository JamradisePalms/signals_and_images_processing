#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <complex>
#include <time.h>
#include <map>

using namespace std;
const double pi = 3.141592;
const double m_pi = 3.141592;

vector<double> sin_signal(double F, double T, double phi, int N) {
    vector<double> signal;
    double dt = 1 / F;
    for (int i = 0; i < F * N; i++) {
        double t = i * dt;
        signal.push_back(sin(2 * pi * t / T + F * phi));
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

vector<double> square_signal(double F, double T, double phi, double t, int N) {
    vector<double> signal(N * F, 0.0);
    double dt = 1.0 / F;
    for (int i = 0; i < F * N; ++i) {
        double t = i * dt;
        signal[i] = sin(2 * pi * t / T + F * phi) > 0? 1.0 : 0;
    }

    return signal;
}

vector<double> double_square_signal(double t, double F, int N) {
    vector<double> signal(N * F, 0.0);
    int pulse1_start = rand() % static_cast<int>(N * F - t * F);
    int pulse2_start;
    do {
        pulse2_start = rand() % static_cast<int>(N * F - t * F);
    } while (abs(pulse2_start - pulse1_start) < t * F);
    for (int i = 0; i < t * F + 1; i++) {
        if (pulse2_start + i >= F * N - 1) return signal;
        signal[pulse1_start + i] = 1.0;
        signal[pulse2_start + i] = 1.0;
    }
    return signal;
}

void save_to_csv(const vector<uint16_t>& signal, double F, int N, const string& filename) {
    double q_min = -1.0;
    double q_max = 1.0;
    int q_levels = 1024;
    
    std::ofstream csv_file;
    csv_file.open(filename);
    csv_file << "time,signal\n";
    for (int i = 0; i < F * N; i++) {
        double signal_val = (q_max - q_min) * double(signal[i]) / (q_levels - 1) + q_min;
        csv_file << (i / F) << "," << signal_val << "\n";
    }
    csv_file.close();
}

vector<uint16_t> quantize(const vector<double>& signal) {
    vector<uint16_t> q_signal(signal.size());
    double q_min = -1.0;
    double q_max = 1.0;
    int q_levels = 1024;
    double step = (q_max - q_min) / (q_levels - 1);
    for (long i = 0; i < signal.size(); i++) {
        double value = max(q_min, min(signal[i], q_max)); // [q_min, q_max]
        q_signal[i] = static_cast<uint16_t>((value - q_min) / step);
    }
    return q_signal;
}

void ComplexBitReverse(complex<double>* data, int size) {
    int middle = size/2, revSize = size - 1, j = 0;
    for (int i = 0; i < revSize; ++i) {
        if (i < j) swap(data[i], data[j]);
        int k = middle;
        while (k <= j) {
            j -= k;
            k /= 2;
        }
        j += k;
    }
}

void FftDit(complex<double>* data, int size, int sizeLog2, int dir ) {
    ComplexBitReverse(data, size);
    int ptsInLeftDft,ptsInRightDft = 1;
    for (int stage = 1; stage <= sizeLog2; ++stage) {
        ptsInLeftDft = ptsInRightDft;
        ptsInRightDft *= 2;
        complex<double> twiddle = complex<double>(1.0, 0.0);
        double trigArg = m_pi / ptsInLeftDft;
        complex<double> wFactor = complex<double>(cos(trigArg),-sin(trigArg)*dir);

        for(int butterflyPos = 0; butterflyPos < ptsInLeftDft; ++butterflyPos) {
            for(int topNode=butterflyPos; topNode < size; topNode+=ptsInRightDft) {
                int botNode = topNode + ptsInLeftDft;
                complex<double> temp = data[botNode] * twiddle;
                data[botNode] = data[topNode] - temp;
                data[topNode] += temp;
            }
            twiddle *= wFactor;
        }
    }
}

void save_to_csv_magnitude(const vector<double>& signal, double F, int N, const string& filename) {
    std::ofstream csv_file;
    csv_file.open(filename);
    csv_file << "magnitude\n";
    for (int i = 0; i < F * N; i++) {
        csv_file << signal[i] << "\n";
    }
    csv_file.close();
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

bool is_power_of_two(int n) {
    return ((n & (n - 1)) == 0);
}

int next_power_of_two(int n) {
    if (n <= 1) return 1;
    return pow(2, ceil(log2(n)));
}

vector<double> pad_to_power_of_two(const vector<double>& signal) {
    int original_size = signal.size();

    if (is_power_of_two(original_size)) {
        return signal;
    }

    int new_size = next_power_of_two(original_size);
    vector<double> padded_signal(new_size, 0.0);
    
    for (int i = 0; i < original_size; ++i) {
        padded_signal[i] = signal[i];
    }
    
    return padded_signal;
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

double box_muller() {
    double s;
    double x, y;
    do {
        x = static_cast<double>(rand()) / RAND_MAX * 2 - 1;
        y = static_cast<double>(rand()) / RAND_MAX * 2 - 1;
        s = x * x + y * y;
    } while (abs(s) < 1e-10 || s > 1);
    double z0 = x * sqrt(-2 * log(s) / s);
    return z0;
}

vector<double> gaussian_noise(const vector<double>& signal, double S, double M) {
    vector<double> noisy_signal = signal;
    for (int i = 0; i < signal.size(); i++) {
        noisy_signal[i] += box_muller() * sqrt(S) + M;
    }
    return noisy_signal;
}

void histogram(const vector<double>& signal, int num_bins, const string& filename) {
    map<int, int> hist;

    double min_value = *min_element(signal.begin(), signal.end());
    double max_value = *max_element(signal.begin(), signal.end());
    double bin_width = (max_value - min_value) / num_bins;

    for (double value : signal) {
        int bin_index = static_cast<int>((value - min_value) / bin_width);
        hist[bin_index] += 1;
    }
    
    std::ofstream csv_file;
    csv_file.open(filename);
    csv_file << "bin" << "," << "count\n";
    for (int i = 0; i < num_bins; i++) {
        csv_file << i << "," << hist[i] << "\n";
    }
}


int main()
{
    srand(time(0));
    double T = 2.0;
    double phi = 0.0;
    double F = 48.0;
    double t = 1;
    int N = 6;
    double P = 0.05;
    double S = 0.01;
    double M = 0.0;
    int num_bins = 50;
    // Process sine signal
    vector<double> sin_sig = sin_signal(F, T, phi, N);
    sin_sig = pad_to_power_of_two(sin_sig);
    vector<uint16_t> q_sin_sig = quantize(sin_sig);
    save_to_csv(q_sin_sig, F, N, "sin_signal_sampled.csv");
    // vector<complex<double>> sin_sig_complex(sin_sig.size());
    // transform(sin_sig.begin(), sin_sig.end(), sin_sig_complex.begin(), [](double x) { return complex<double>(x, 0.0); });
    // FftDit(sin_sig_complex.data(), sin_sig_complex.size(), log2(sin_sig_complex.size()), 1);
    // vector<double> sin_magnitude(sin_sig_complex.size());
    // transform(sin_sig_complex.begin(), sin_sig_complex.end(), sin_magnitude.begin(), [](complex<double> x) { return abs(x); });
    // save_to_csv_magnitude(sin_magnitude, F, N, "sin_signal_complex.csv");

    // Process triangle signal
    vector<double> triangle_sig = triangle_signal(F, T, phi, t, N);
    triangle_sig = pad_to_power_of_two(triangle_sig);
    vector<uint16_t> q_triangle_sig = quantize(triangle_sig);
    save_to_csv(q_triangle_sig, F, N, "triangle_signal_sampled.csv");
    // vector<complex<double>> triangle_sig_complex(triangle_sig.size());
    // transform(triangle_sig.begin(), triangle_sig.end(), triangle_sig_complex.begin(), [](double x) { return complex<double>(x, 0.0); });
    // FftDit(triangle_sig_complex.data(), triangle_sig_complex.size(), log2(triangle_sig_complex.size()), 1);
    // vector<double> triangle_magnitude(triangle_sig_complex.size());
    // transform(triangle_sig_complex.begin(), triangle_sig_complex.end(), triangle_magnitude.begin(), [](complex<double> x) { return abs(x); });
    // save_to_csv_magnitude(triangle_magnitude, F, N, "triangle_signal_complex.csv");

    // Process square signal
    vector<double> square_sig = square_signal(F, T, phi, t, N);
    square_sig = pad_to_power_of_two(square_sig);
    vector<uint16_t> q_square_sig = quantize(square_sig);
    save_to_csv(q_square_sig, F, N, "square_signal_sampled.csv");
    // vector<complex<double>> square_sig_complex(square_sig.size());
    // transform(square_sig.begin(), square_sig.end(), square_sig_complex.begin(), [](double x) { return complex<double>(x, 0.0); });
    // FftDit(square_sig_complex.data(), square_sig_complex.size(), log2(square_sig_complex.size()), 1);
    // vector<double> square_magnitude(square_sig_complex.size());
    // transform(square_sig_complex.begin(), square_sig_complex.end(), square_magnitude.begin(), [](complex<double> x) { return abs(x); });
    // save_to_csv_magnitude(square_magnitude, F, N, "square_signal_complex.csv");

    // Process double square signal
    vector<double> double_square_sig = double_square_signal(t, F, N);
    double_square_sig = pad_to_power_of_two(double_square_sig);
    vector<uint16_t> q_double_square_sig = quantize(double_square_sig);
    save_to_csv(q_double_square_sig, F, N, "double_square_signal_sampled.csv");
    vector<complex<double>> double_square_sig_complex(double_square_sig.size());
    transform(double_square_sig.begin(), double_square_sig.end(), double_square_sig_complex.begin(), [](double x) { return complex<double>(x, 0.0); });
    FftDit(double_square_sig_complex.data(), double_square_sig_complex.size(), log2(double_square_sig_complex.size()), 1);
    vector<double> double_square_magnitude(double_square_sig_complex.size());
    transform(double_square_sig_complex.begin(), double_square_sig_complex.end(), double_square_magnitude.begin(), [](complex<double> x) { return abs(x); });
    save_to_csv_magnitude(double_square_magnitude, F, N, "double_square_signal_complex.csv");

    //Creating noisy impulse signal
    vector<double> noisy_sin_sig = impulse_noise(sin_sig, P);
    save_to_csv_double(noisy_sin_sig, F, N, "noisy_sin_signal.csv");

    //Creating noisy triangle signal
    vector<double> noisy_triangle_sig = impulse_noise(triangle_sig, P);
    save_to_csv_double(noisy_triangle_sig, F, N, "noisy_triangle_signal.csv");

    //Creating noisy square signal
    vector<double> noisy_square_sig = impulse_noise(square_sig, P);
    save_to_csv_double(noisy_square_sig, F, N, "noisy_square_signal.csv");

    //Creating noisy double square signal
    vector<double> noisy_double_square_sig = impulse_noise(double_square_sig, P);
    save_to_csv_double(noisy_double_square_sig, F, N, "noisy_double_square_signal.csv");

    // GAUSS
    //Creating noisy gaussian signal
    vector<double> noisy_gaussian_sin_sig = gaussian_noise(sin_sig, S, M);
    save_to_csv_double(noisy_gaussian_sin_sig, F, N, "noisy_gaussian_sin_signal.csv");

    //Creating noisy gaussian triangle signal
    vector<double> noisy_gaussian_triangle_sig = gaussian_noise(triangle_sig, S, M);
    save_to_csv_double(noisy_gaussian_triangle_sig, F, N, "noisy_gaussian_triangle_signal.csv");

    //Creating noisy gaussian square signal
    vector<double> noisy_gaussian_square_sig = gaussian_noise(square_sig, S, M);
    save_to_csv_double(noisy_gaussian_square_sig, F, N, "noisy_gaussian_square_signal.csv");

    //Creating noisy gaussian double square signal
    vector<double> noisy_gaussian_double_square_sig = gaussian_noise(double_square_sig, S, M);
    save_to_csv_double(noisy_gaussian_double_square_sig, F, N, "noisy_gaussian_double_square_signal.csv");

    //Histograms
    histogram(noisy_gaussian_sin_sig, num_bins, "histogram_noisy_gaussian_sin_signal.csv");
    histogram(noisy_gaussian_triangle_sig, num_bins, "histogram_noisy_gaussian_triangle_signal.csv");
    histogram(noisy_gaussian_square_sig, num_bins, "histogram_noisy_gaussian_square_signal.csv");
    histogram(noisy_gaussian_double_square_sig, num_bins, "histogram_noisy_gaussian_double_square_signal.csv");

    return 0;
}