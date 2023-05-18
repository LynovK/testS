//
// Created by perseverance on 09.04.2023.
//

#ifndef SLAE_FIVE_DIAG_H
#define SLAE_FIVE_DIAG_H

#include <iostream>
#include <vector>

/*template<typename T>
struct element {
    unsigned int i;
    unsigned int j;
    T v;
};*/

template<typename T>
class five_diag {
    const T a, b;
    const unsigned int L;
public:
    std::vector<element<double>> matrix;

    five_diag(const unsigned int n, const T a, const T b, const unsigned int L) noexcept: a(a), b(b), L(L) {
        std::vector<element<double>> matrix_(5 * n - 2 * L - 2);
        for (unsigned int i = 0; i < n; ++i) {
            matrix_.push_back({i, i, 2 * b});
        }
        for (unsigned int i = 0; i < n - 1; ++i) {
            matrix_.push_back({i + 1, i, a});
            matrix_.push_back({i, i + 1, a});
        }
        for (unsigned int i = 0; i < n - L; ++i) {
            matrix_.push_back({i + L, i, a});
            matrix_.push_back({i, i + L, a});
        }
        matrix = matrix_;
    }

    T lambda_max() const noexcept {
        return 2 * (b + 2 * a * cos(M_PI / (L + 1)));
    }

    T lambda_min() const noexcept {
        return 2 * (b - 2 * a * cos(M_PI / (L + 1)));
    }
};

#endif //SLAE_FIVE_DIAG_H