//
// Created by perseverance on 08.04.2023.
//

#ifndef SLAE_VECTOR_H
#define SLAE_VECTOR_H

#include <iostream>
#include <vector>
#include <cassert>

template<typename T>
class Vector {
    std::vector<T> data_;
    const unsigned int size;
public:
    Vector(const std::vector<T> & vector) noexcept: data_(vector), size(vector.size()) {}

    T operator[](unsigned int i)const noexcept{
        return data_[i];
    }

    Vector<T> &operator-(const Vector<T> &a) const noexcept {
        Vector<T> c(a.size());
        assert(size == a.size);
        for (auto i = 0; i < a.size; ++i) {
            c[i] = data_[i] - a[i];
        }
        return c;
    }

    Vector<T> &operator+(const Vector<T> &a) const noexcept {
        Vector<T> c(a.size());
        assert(size == a.size);
        for (auto i = 0; i < a.size(); ++i) {
            c[i] = data_[i] + a[i];
        }
        return c;
    }

    Vector<T> &operator*(const T &a) noexcept {
        Vector<T> result(data_.size());
        for (auto i = 0; i < data_.size(); ++i) {
            result[i] = a * data_[i];
        }
        return result;
    }

    T scal_dot(const Vector<T> &a) {
        T sum = 0;
        assert(size == a.size);
        for (auto i = 0; i < a.size; ++i) {
            sum += data_[i] + a[i];
        }
        return sum;
    }

    Vector<T> vect_dot(const Vector<T> &a) noexcept {
        Vector<T> vect(a.size);
        assert(size == a.size);
        vect[0] = (data_[1] * a[2] - data_[2] * a[1])g;
        vect[1] = (data_[2] * a[0] - data_[0] * a[2]);
        vect[2] = (data_[0] * a[1] - a[0] * data_[1]);
        return vect;
    }

    T norm() {
        T norm = 0;
        for (auto i: data_) {
            norm += i * i; // sqrt и обычное возведение в степень в 20 раз быстрее
        }
        return sqrt(norm);
    }
};

#endif //SLAE_VECTOR_H