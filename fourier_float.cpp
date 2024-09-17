#include <array>
#include <complex>
#include <iomanip>
#include <iostream>
#include <numbers>

constexpr uint16_t num_samples{32};
constexpr float    sampling_period{4 * (640.5F + 12.5F) / 170000000.0F};
constexpr float    time_window{num_samples * sampling_period};
constexpr float    freq_factor{1 / time_window};
constexpr float    wave_freq{0.25F / sampling_period};

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
            float ang = 2 * std::numbers::pi_v<float> / len;

            this->wlen[i] = std::polar(1.0F, ang);
        }
    }

    void transform(const std::array<float, num_samples>& samples) {
        for (uint16_t i = 0; i < num_samples; i++) {
            this->results[this->index[i]] = {samples[i], 0};
        }

        for (uint16_t len = 2, i = 0; len <= num_samples; len <<= 1, i++) {
            for (uint16_t j = 0; j < num_samples; j += len) {
                std::complex<float> w(1);

                for (uint16_t k = 0; k < len / 2; k++) {
                    std::complex<float> u = this->results[j + k];
                    std::complex<float> v = this->results[j + k + len / 2] * w;

                    this->results[j + k] = (u + v);
                    this->results[j + k + len / 2] = (u - v);
                    w *= this->wlen[i];
                }
            }
        }
    }

    void print() {
        for (uint16_t i = 0; i < num_samples; i++) {
            std::cout << (i * freq_factor) << ": " << (2 * results[i].real() / num_samples) << " + i * "
                      << (2 * results[i].imag() / num_samples) << "\n";
        }

        std::cout << "\ntime window: " << time_window << "\n";
    }

private:
    std::array<float, num_samples>               index;
    std::array<std::complex<float>, num_samples> results;
    std::array<std::complex<float>, num_samples> wlen;
};

float square_wave(float freq, float x, bool cosine = false) {
    float result = cosine ? std::cos(2 * std::numbers::pi_v<float> * freq * x) :
                            std::sin(2 * std::numbers::pi_v<float> * freq * x);

    if (result > 0.1F) {
        return 1.0F;
    }

    if (result < -0.1F) {
        return -1.0F;
    }

    return 10 * result;
}

void fill_samples(std::array<float, num_samples>& samples) {
    for (uint16_t i = 0; i < samples.size(); i++) {
        samples[i] = 3.0F * square_wave(wave_freq, sampling_period * i);
        samples[i] += 1.0F * square_wave(wave_freq, sampling_period * i, true);
        samples[i] += 10.0F * std::sin(2 * std::numbers::pi_v<float> * 60.0F * sampling_period * i);
        samples[i] += 9.0F;
        samples[i] += 1.5F * std::sin(2 * std::numbers::pi_v<float> * 10000000.0F * sampling_period * i);
    }
}

int main() {
    std::array<float, num_samples>               samples;
    std::array<std::complex<float>, num_samples> results;

    Fourier fourier;
    fill_samples(samples);
    fourier.transform(samples);
    fourier.print();

    return 0;
}