//
// Created by perseverance on 20.02.23.
//
#include <gtest/gtest.h>
#include "../src/Tridiagonal/ThreeDiagonalSolver.h"


TEST(ThreeDiagTest, firsttest) {
    std::vector<elements<double>> example = {{0,2,1}, {1,4,2}, {5,8,1}, {3,9,1}, {7,2,0}};
    std::vector<double> d = {3, 7, 14, 13, 9};
    std::vector<double> x = {1, 1, 1, 1, 1, 1};

    ThreeDiagonalMatrix<double> first(example);


    std::vector<double> solve;
    solve = solver(first, d);

    for (int i = 0; i < solve.size(); i++){
        ASSERT_NEAR(solve[i], x[i], 0.001);
    }
}


TEST(ThreeDiagTest, secondtest) {
    std::vector<elements<double>> example = {{0,2,1}, {2,7,1}, {3,9,1}, {5,6,0}};
    std::vector<double> d = {4, 20, 43, 26};
    std::vector<double> x = {1, 2, 4, 1};

    ThreeDiagonalMatrix<double> first(example);


    std::vector<double> solve;
    solve = solver(first, d);

    for (int i = 0; i < solve.size(); i++){
        ASSERT_NEAR(solve[i], x[i], 0.000001);
    }
}

