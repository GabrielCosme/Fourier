#include <array>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>

constexpr uint16_t num_samples{64};
constexpr double   sampling_period{653.0 / 170000000.0};
constexpr double   time_window{num_samples * sampling_period};
constexpr double   freq_factor{1 / time_window};
constexpr double   wave_freq{0.25 / sampling_period};

class Fourier {
public:
    Fourier() {
        for (uint16_t i = 0; i < num_samples; i++) {
            this->index[i] = i;
        }

        for (uint16_t i = 1, j = 0; i < num_samples; i++) {
            uint16_t bit = num_samples >> 1;

            for (; j & bit; bit >>= 1) {
                j ^= bit;
            }

            j ^= bit;

            if (i < j) {
                std::swap(this->index[i], this->index[j]);
            }
        }

        for (uint16_t len = 2, i = 0; len <= num_samples; len <<= 1, i++) {
            double ang = 2 * std::numbers::pi / len;

            this->wlen[i] = std::polar(1.0, ang);
        }
    }

    void transform(const std::array<double, num_samples>& samples) {
        for (uint16_t i = 0; i < num_samples; i++) {
            this->results[this->index[i]] = {samples[i], 0};
        }

        for (uint16_t len = 2, i = 0; len <= num_samples; len <<= 1, i++) {
            for (uint16_t j = 0; j < num_samples; j += len) {
                std::complex<double> w(1);

                for (uint16_t k = 0; k < len / 2; k++) {
                    std::complex<double> u = this->results[j + k];
                    std::complex<double> v = this->results[j + k + len / 2] * w;

                    this->results[j + k] = (u + v);
                    this->results[j + k + len / 2] = (u - v);
                    w *= this->wlen[i];
                }
            }
        }
    }

    void print() {
        for (uint16_t i = 0; i < num_samples; i++) {
            std::cout << (i * freq_factor) << ": " << (2 * std::abs(results[i]) / num_samples) << "\n";
        }

        std::cout << "\ntime window: " << time_window << "\n";
    }

private:
    std::array<double, num_samples>               index;
    std::array<std::complex<double>, num_samples> results;
    std::array<std::complex<double>, num_samples> wlen;
};

void fill_samples(std::array<double, num_samples>& samples) {
    for (uint16_t i = 0; i < samples.size(); i++) {
        samples[i] = 3 * sin(2 * std::numbers::pi * wave_freq * sampling_period * i);
        samples[i] += 10 * sin(2 * std::numbers::pi * 60 * sampling_period * i);
        samples[i] += 9;
        samples[i] += 1.5 * sin(2 * std::numbers::pi * 10000000 * sampling_period * i);
    }
}

int main() {
    std::array<double, num_samples>               samples;
    std::array<std::complex<double>, num_samples> results;

    Fourier fourier;
    fill_samples(samples);
    fourier.transform(samples);
    fourier.print();

    return 0;
}