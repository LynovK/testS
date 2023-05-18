#ifndef SLAE_CHEBISHEVMPI_H
#define SLAE_CHEBISHEVMPI_H

#include <iostream>
#include "../Iteration_methods/Iteration_methods.h"
#include <cmath>


template<typename T>
std::vector<T> MPI_chebishev(const CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b,
                             const std::vector<T> &x_0, const std::vector<T> &tau) noexcept {
    unsigned int count = 0;
    std::vector<T> x = x_0;
    std::ofstream out;
    out.open("5.txt");

    std::vector<T> r = CSR.dot(x) - b; //Нач приближение и невязка
    while (tolerance < r_inf<T>(CSR, x, b)) {
        count++;
        for (auto i = 0; i < tau.size(); ++i) {
            x = x - r * tau[i];
            r = CSR.dot(x) - b;
            out << r_inf(CSR, x, b) << ' ';
        }
    }
    std::cout << "Chebishev count: " << count << std::endl;
    return x;
}

template<typename T>
std::vector<T> Chebishev_solver(unsigned int n, unsigned int r, T lambda_min, T lambda_max) {
    std::vector<T> z;
    z.resize(n);
    std::vector<T> betta;
    betta.resize(n);
    const double pi = M_PI;

    //Корни полинома Чебышева
    betta[0] = sin(pi / (2 * n));
    T sin_alpha = sin(pi / n);
    T cos_alpha = cos(pi / n);
    z[0] = cos(pi / (2 * n));
    //z[1] = z[0] * cos_alpha - betta[0] * sin_alpha;

    for (auto i = 0; i < n - 1; ++i) {
        betta[i + 1] = betta[i] * cos_alpha + z[i] * sin_alpha;
        z[i + 1] = z[i] * cos_alpha - betta[i] * sin_alpha;
    }

    T c = (lambda_min + lambda_max) / 2;
    T k = (lambda_max - lambda_min) / 2;
    //Афинное преобразование
    for (auto i = 0; i < z.size(); ++i) {
        z[i] = c + k * z[i];
    }

    //Переупорядоченный массив индексов
    std::vector<unsigned int> indexes = {0, 1};
    std::vector<unsigned int> ind;
    ind.resize(n);
    indexes.resize(n);
    unsigned int point = 4;
    for (auto i = 1; i < r; ++i) {
        for (auto j = 0; j < point / 2; ++j) {
            ind[2 * j] = indexes[j];
            ind[2 * j + 1] = point - 1 - indexes[j];
            //indexes = ind;
        }
        indexes = ind;
        point *= 2;
    }

    //Вектор тау
    std::vector<T> tau;
    tau.resize(n);
    for (auto i = 0; i < z.size(); ++i) {
        tau[i] = 1 / z[indexes[i]];
    }
    return tau;
}

#endif //SLAE_CHEBISHEVMPI_H
