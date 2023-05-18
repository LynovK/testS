//
// Created by perseverance on 29.04.2023.
//

#ifndef SLAE_DENSE_MATRIX_H
#define SLAE_DENSE_MATRIX_H

#include <iostream>
#include <vector>

template<typename T>
class dense_matrix {
private:
    std::vector<T> matrix;
    unsigned int n;
    unsigned int m;
public:
    dense_matrix(const std::vector<T> &vect, const unsigned int M, const unsigned int N) noexcept: matrix(vect), n(N),
                                                                                                   m(M) {}

    std::vector<T> operator*(const std::vector<T> x) const noexcept {
        assert(x.size() == n);
        std::vector<T> res;
        res.reserve(m);
        for (size_t i = 0; i < m; ++i) {
            T sum = 0;
            for (auto j = 0; j < n; ++j) {
                sum += matrix[m * i + j] * x[j];
            }
            res.push_back(sum);
        }
        return res;
    }

    T& operator()(unsigned int i, unsigned int j) {
        assert(i < m && j < n);
        return matrix[3 * i + j];
    }

    void operator()(unsigned int i, unsigned int j, T k) {
        assert(i < m && j < n);
        matrix[3 * i + j] = k;
    }
};

#endif //SLAE_DENSE_MATRIX_H