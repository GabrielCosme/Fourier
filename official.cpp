#include <complex>
#include <numbers>
#include <array>

void fft(const std::array<double, 1024>& samples, std::array<std::complex<double>, 1024>& results) {
    int n = samples.size();

    for (int i = 0; i < n; i++) {
        results[i] = {samples[i], 0};
    }

    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;

        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }

        j ^= bit;

        if (i < j) {
            std::swap(results[i], results[j]);
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * std::numbers::pi / len;

        std::complex<double> wlen(cos(ang), sin(ang));

        for (int i = 0; i < n; i += len) {
            std::complex<double> w(1);

            for (int j = 0; j < len / 2; j++) {
                std::complex<double> u = results[i + j];
                std::complex<double> v = results[i + j + len / 2] * w;

                results[i + j] = (u + v);
                results[i + j + len / 2] = (u - v);
                w *= wlen;
            }
        }
    }
}