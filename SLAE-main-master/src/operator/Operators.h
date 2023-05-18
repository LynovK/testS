//
// Created by perseverance on 08.04.2023.
//

#ifndef SLAE_OPERATORS_H
#define SLAE_OPERATORS_H

#include <iostream>
#include <vector>
#include <cassert>

template<typename T>
struct basis {
    std::vector<T> I;// = {1, 0, 0};
    std::vector<T> J;// = {0, 1, 0};
    std::vector<T> K;// = {0, 0, 1};
};

template<typename T>
std::vector<T> operator+(const std::vector<T> &a, const std::vector<T> &b) {
    std::vector<T> sum(a.size());
    assert(a.size() == b.size());
    for (auto i = 0; i < a.size(); ++i) {
        sum[i] += a[i] + b[i];
    }
    return sum;
}

template<typename Type>
std::vector<Type> operator-(const std::vector<Type> &a, const std::vector<Type> &b) {
    std::vector<Type> result(a.size());
    assert(a.size() == b.size());
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
T scal_dot(const std::vector<T> &a, const std::vector<T> &b) noexcept {
    T scal = 0;
    assert(a.size() == b.size());
    for (auto i = 0; i < a.size(); ++i) {
        scal += a[i] * b[i];
    }
    return scal;
}

template<typename T>
std::vector<T> vect_dot(const std::vector<T> &a, const std::vector<T> &b) noexcept {
    std::vector<T> vect(a.size());
    vect[0] = (a[1] * b[2] - a[2] * b[1]);
    vect[1] = (a[2] * b[0] - a[0] * b[2]);
    vect[2] = (a[0] * b[1] - b[0] * a[1]);
    return vect;
}

template<typename T>
std::vector<T> vect_dot(const std::vector<T> &a, const std::vector<T> &b, const basis<T> &basis) noexcept {
    std::vector<T> vect(a.size());
    vect = (a[1] * b[2] - a[2] * b[1]) * basis.I + (a[2] * b[0] - a[0] * b[2]) * basis.J +
           (a[0] * b[1] - b[0] * a[1]) * basis.K;
    return vect;
}

#endif //SLAE_OPERATORS_H