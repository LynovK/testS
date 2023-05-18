//
// Created by perseverance on 13.02.23.
//
#pragma once
#ifndef SLAE_COMPRESSEDMATRIX_H
#define SLAE_COMPRESSEDMATRIX_H

#include <iostream>
#include <vector>


template<typename T>
struct element {
    unsigned int i;
    unsigned int j;
    T v;
};

template<typename T>
class CompressedMatrix {
private:
    std::vector<T> vector_of_values_;
    std::vector<unsigned int> colomn_index_vector_;
    std::vector<unsigned int> counts_vector_;
public:
    CompressedMatrix<T>(std::vector<element<T>> matrix_, unsigned int n,
                        unsigned int m) noexcept { // matrix приходит уже отсортированная
        std::vector<T> vector_of_values;
        std::vector<unsigned int> colomn_index_vector;
        std::vector<unsigned int> counts_vector;

        counts_vector.resize(n + 1);
        counts_vector[0] = 0;

        for (auto i = 0; i < matrix_.size(); ++i) {
            vector_of_values.push_back(matrix_[i].v);         //Получаем строку с ненулевыми элементами матрицы
            colomn_index_vector.push_back(matrix_[i].j);      //Получаем строку с номерами столбцов
            counts_vector[matrix_[i].i + 1] += 1;
        }

        for (auto i = 1; i < counts_vector.size(); ++i) {
            counts_vector[i] += counts_vector[i - 1];

        }

        vector_of_values_ = vector_of_values;
        colomn_index_vector_ = colomn_index_vector;
        counts_vector_ = counts_vector;
    }

    const std::vector<T> &getValue() const noexcept {
        return vector_of_values_;
    }

    const std::vector<unsigned int> &get_index() const noexcept {
        return colomn_index_vector_;
    }

    const std::vector<unsigned int> &get_count() const noexcept {
        return counts_vector_;
    }

    //Умножение на вектор
    std::vector<T> dot(const std::vector<T> &X) const noexcept {
        std::vector<T> x(X.size(), 0);
        for (auto i = 0; i < counts_vector_.size() - 1; ++i) {
            for (auto j = counts_vector_[i]; j < counts_vector_[i + 1]; ++j) {
                x[i] += vector_of_values_[j] * X[colomn_index_vector_[j]];
            }
        }
        return x;
    }

    //взятие по индексу
    T operator()(unsigned int i, unsigned int j) const noexcept {
        for (auto k = counts_vector_[i]; k < counts_vector_[i + 1]; ++k) {
            if (colomn_index_vector_[k] == j)
                return vector_of_values_[k];
        }
        return 0;
    }
};

#endif //SLAE_COMPRESSEDMATRIX_H