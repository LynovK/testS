#ifndef SLAE_RELATED_DIRECTIONS_METHOD_H
#define SLAE_RELATED_DIRECTIONS_METHOD_H

#include "Iteration_methods.h"
#include <ctime>

template<typename T>
T scal_dot(const std::vector<T> &a, const std::vector<T> &b) noexcept {
    T scal = 0;
    if (a.size() == b.size()) {
        for (auto i = 0; i < a.size(); ++i) {
            scal += a[i] * b[i];
        }
    } else { std::cout << "Different size of vectors in scal dot" << std::endl; }
    return scal;
}

template<typename T>
std::vector<T> Related_directions(const CompressedMatrix<T> &CSR, const T tolerance, const std::vector<T> &b,
                                  const std::vector<T> &x_0) noexcept {
    unsigned int count = 0;
    std::vector<T> x = x_0;
    std::vector<T> r = CSR.dot(x) - b; //Нач приближение и невязка
    std::vector<T> d = r;
    std::ofstream out;

    //std::vector<T> rk;
    T alfa;
    out.open("3.txt");

    while (tolerance < r_inf<T>(CSR, x, b)) {
        count++;
        alfa = scal_dot(d, r) / scal_dot(d, CSR.dot(d));
        x = x - d * alfa;
        //rk = r;
        const T scal = scal_dot(d, r);
        r = CSR.dot(x) - b;
        d = r + d * (scal_dot(r, r) / scal);
        out << r_inf(CSR, r, b) << ' ';
    }
    std::cout << "Related directions count: " << count << std::endl;
    //std::cout << "runtime: " << clock()/1000.0 << std::endl;
    return x;
}

#endif //SLAE_RELATED_DIRECTIONS_METHOD_H
