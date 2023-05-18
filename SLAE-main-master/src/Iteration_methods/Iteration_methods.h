//
// Created by perseverance on 05.03.23.
//

#pragma once
#ifndef SLAE_ITERATION_METHODS_H
#define SLAE_ITERATION_METHODS_H

#include <cmath>
#include "../Compressed_spars_row/CompressedMatrix.h"
#include "../Compressed_spars_row/CompressedSorting.h"
#include <fstream>

template<typename Type>
std::vector<Type> operator-(const std::vector<Type> &a, const std::vector<Type> &b) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

template<typename Type>
std::vector<Type> operator*(const std::vector<Type> &a, const Type &c) {
    std::vector<Type> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = c * a[i];
    }
    return result;
}

template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> sum(a.size());
    for (auto i = 0; i < a.size(); ++i) {
        sum[i] += a[i] + b[i];
    }
    return sum;
}

/*template<typename T>
T get_r(const CompressedMatrix<T> &CSR, const std::vector<T> &x, const std::vector<T> &b) noexcept {
    T r = 0;
    for (size_t j = 0; j < b.size(); ++j) {
        r += (b[j] - CSR.dot(x)[j]) * (b[j] - CSR.dot(x)[j]); //евклидова длина стобца невязки
    }
    return r;
}*/

template<typename T>
T r_inf(const CompressedMatrix<double> &CSR, const std::vector<T> &x, const std::vector<T> &b) noexcept {
    T r = 0;
    for (size_t j = 0; j < b.size(); ++j) { //size_t = unsigned int
        if (fabs(b[j] - CSR.dot(x)[j]) > r)
            r = fabs(b[j] - CSR.dot(x)[j]);
    }
    return r;
}

template<typename T>
std::vector<T> MPI(CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b,
                   const std::vector<T> &x_0, const T tau) noexcept {
    unsigned int count = 0;
    std::ofstream out;
    out.open("1.txt");
    std::vector<T> x = x_0; //Нач приближение и невязка
    std::vector<T> r = CSR.dot(x) - b;
    while (tolerance < r_inf<T>(CSR, x, b)) {
        count++;
        x = x - r * tau;
        r = CSR.dot(x) - b;
        out << r_inf(CSR, x, b) << ' ';
    }
    std::cout << "MPI count: " << count << std::endl;
    return x;
}

/*template<typename T>
std::vector<T> Jacoby(CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b,
                      const std::vector<T> &x0) {
    std::vector<T> xk = x0, xt(x0.size(), 0);
    unsigned long long count = 0;
    while (tolerance < r_inf<T>(CSR, xk, b)) {
        count++;
        for (size_t k = 0; k < x0.size(); k++) {
            xt[k] = (1 / CSR.element(k, k)) * (b[k] - CSR.dot(xk)[k] + CSR.element(k, k) * xk[k]);
        }
        xk = xt;
    }
    return xk;
}*/

template<typename T>
std::vector<T>
Jacoby(const CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b, const std::vector<T> &x0) {
    std::vector<T> x_0 = x0;
    std::vector<T> x1(x0.size());
    unsigned int count = 0;
    while (tolerance < r_inf<T>(CSR, x1, b)) {
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
            //x1[i] = 1 / diag * (b[i] - sum);
        }
        x1 = x_0;
    }
    std::cout << "Jacoby count: " << count << std::endl;
    return x_0;
}

template<typename T>
std::vector<T>
G_Z(const CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b, const std::vector<T> &x0) {
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
    }
    std::cout << "G_Z count: " << count << std::endl;
    return x_0;
}

template<typename T>
std::vector<T>
SOR(const CompressedMatrix<T> &CSR, const T tolerance, const T omega, const std::vector<T> &b,
    const std::vector<T> &x0) {
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
                    sum += omega * CSR.getValue()[j] * x_0[CSR.get_index()[j]];
                }
            }
            x_0[i] = 1 / diag * (omega * b[i] - sum - (omega - 1) * diag * x_0[i]);
        }
    }
    std::cout << "SOR count: " << count << std::endl;
    return x_0;
}

#endif //SLAE_ITERATION_METHODS_H