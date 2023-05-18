//
// Created by perseverance on 29.04.2023.
//

#ifndef SLAE_GMRES_H
#define SLAE_GMRES_H

#include <iostream>
#include "../operator/Operators.h"
#include "../Matrices/dense_matrix.h"
#include "../Compressed_spars_row/CompressedMatrix.h"
#include "../Compressed_spars_row/CompressedSorting.h"

template<typename T>
struct Krilov_basis {
    std::vector<std::vector<T>> basis;
};

template <typename T>
/*template<typename T>
std::vector<T>
//dense_matrix<T>
GMRES(const CompressedMatrix<T> &CSR, const std::vector<T> &x_0, const std::vector<T> &b, unsigned int n,
      unsigned int m) {
    std::vector<T> x = x_0;
    std::vector<T> r = CSR.dot(x) - b;

    std::vector<T> h((n + 1) * n, 0); //поле миатрицы Гессенберга
    dense_matrix<T> H(h, n + 1, n); //матрица гессенберга

    std::vector<std::vector<T>> V;
    V.reserve(); //////////////////////// сколько?
    std::vector<T> v_0 = r * (1 / norm(r));
    V.push_back(v_0);

    std::vector<T> t;
    std::vector<T> sinn;
    sinn.reserve();
    std::vector<T> coss;
    coss.reserve(); ////////////////// сколько?

    std::vector<T> e1(x.size(), 0); //x.size =? n;
    e1[0] = 1;
    for (size_t j = 0; j < n; ++j) {
        T betta = norm(r);
        V[0] = r * (1 / betta);
        std::vector<T> bt = betta * e1;

        for (size_t i = 1; i < m + 1; ++i) {
            std::vector<T> w = CSR.dot(V[i]);

            for (size_t k = 1; k < i + 1; ++k) {
                H(k, j) = scal_dot(V[k], w);
                w = w - H(k, j) * V[k];
            }
            H(i + 1, i) = norm(w);
            V[i + 1] = w * (1 / H(i + 1, i));
            r[1] = H(1, i);

            for (size_t k = 2; k < i + 1; ++k) {
                std::vector<T> gamma = coss(k - 1) * r[k - 1] + sinn(k - 1) * H(k, i);
                r[k] = -sinn(k - 1) * r[k - 1] + coss(k - 1) * H(k, i);
                r[k - 1] = gamma;
            }
            T delta = sqrt(r[i] * r[i] + H(i + 1, i) * H(i + 1, i));
            coss[i] = r[i] / delta;
            sinn[i] = H(i + 1, i) / delta;
            r[i] = coss[i] * r[i] + sinn[i] * H(i + 1, i);
            bt[i + 1] = -sinn[i]*bt[i];
            bt[i] = coss[i]*bt[i];
            double rho = abs(bt[i + 1]);
        }
    }
}
*/
#endif //SLAE_GMRES_H