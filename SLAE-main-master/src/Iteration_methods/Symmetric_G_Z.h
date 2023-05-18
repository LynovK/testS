//
// Created by perseverance on 06.04.2023.
//

#ifndef SLAE_SYMMETRIC_G_Z_H
#define SLAE_SYMMETRIC_G_Z_H

#include "Iteration_methods.h"

template<typename T>
std::vector<T>
Symmentric_G_Z(const CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b, const std::vector<T> &x0) {
    std::vector<T> x_0 = x0;
    unsigned int count = 0;
    while (tolerance < r_inf<T>(CSR, x_0, b)) {
        count++;
        for (auto i = 0; i < CSR.get_count().size() - 1; ++i) {
            T sum = 0;
            T diag;
            for (auto j = CSR.get_count()[i]; j < CSR.get_count()[i + 1]; ++j) {
                if (i == CSR.get_index()[j]) {
                    diag = CSR.getValue()[j];
                } else {
                    sum += CSR.getValue()[j] * x_0[CSR.get_index()[j]];
                }
            }
            x_0[i] = 1 / diag * (b[i] - sum);
        }
        for (auto i =CSR.get_count().size() - 1; i >  0; --i) {
            T sum = 0;
            T diag;
            for (auto j = CSR.get_count()[i]; j < CSR.get_count()[i + 1]; ++j) {
                if (i == CSR.get_index()[j]) {
                    diag = CSR.getValue()[j];
                } else {
                    sum += CSR.getValue()[j] * x_0[CSR.get_index()[j]];
                }
            }
            x_0[i] = 1 / diag * (b[i] - sum);
        }
    }
    std::cout << "Symmetric_G_Z count: " << count << std::endl;
    return x_0;
}

#endif //SLAE_SYMMETRIC_G_Z_H