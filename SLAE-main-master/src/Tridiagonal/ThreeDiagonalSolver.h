//
// Created by perseverance on 20.02.23.
//

#ifndef SLAE_THREEDIAGONALSOLVER_H
#define SLAE_THREEDIAGONALSOLVER_H

#include "ThreeDiagonalMatrix.h"
#include <vector>

template<typename T>
std::vector<T> solver(const ThreeDiagonalMatrix<T> &TDM, const std::vector<T> &d) {
    std::vector<T> p;
    auto N = TDM.size();
    p.resize(N);
    p[0] = -TDM.getString(0).c / TDM.getString(0).b;
    //p.push_back(-TDM.getString(0).c / TDM.getString(0).b);

    std::vector<T> q;
    q.resize(N);
    q[0] = d[0] / TDM.getString(0).b;

    std::vector<T> x;
    x.resize(N);

    for (auto i = 0; i < N - 1; ++i){
        auto znam = TDM.getString(i + 1).a*p[i] + TDM.getString(i + 1).b;
        p[i + 1] = -TDM.getString(i + 1).c/znam;
        q[i + 1] = (d[i + 1] - TDM.getString(i + 1).a * q[i])/znam;
    }

    x.back() = (d[N - 1] - TDM.getString(N - 1).a * q[N - 2])/(TDM.getString(N - 1).a * p[N - 2] + TDM.getString(N - 1).b);

    for(auto i = N - 1; i-- > 0; i){
        x[i] = p[i]*x[i + 1] + q[i];
    }

    return x;
}

#endif //SLAE_THREEDIAGONALSOLVER_H