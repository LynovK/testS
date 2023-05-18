#pragma once
#ifndef SLAE_COMPRESSEDSORTING_H
#define SLAE_COMPRESSEDSORTING_H
#include "CompressedMatrix.h"
#include <iostream>
#include <algorithm>
#include <vector>

template<typename T>
void sort_me_plz(std::vector<element<T>> &matrix_){
    std::sort(matrix_.begin(), matrix_.end(), [](const element<T> &first_element, const element<T> &second_element ){
        if (first_element.i != second_element.i) return first_element.i < second_element.i;
        if (first_element.i == second_element.i) return first_element.j < second_element.j;
    });
}

#endif //SLAE_COMPRESSEDSORTING_H