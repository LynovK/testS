#include "../src/Iteration_methods/Grad_spusk.h"
#include <gtest/gtest.h>

TEST(Gradient, Grad_spusk_test) {

    std::vector<element<double>> matrix_CSR = {
            {3, 3, 0},
            {2, 5, 2},
            {1, 4, 1},
            {3, 5, 2},
            {2, 1, 10},
            {12, 2, 3},
            {0, 0, 1},
            {6, 1, 3},
            {5, 1, 10}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res = CompressedMatrix(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 1, 1};
    std::vector<double> x0 = {0, 0, 0};

    double tau = 0.0001;
    double tolerance1 = 1e-10;
    double tolerance2 = 1e-13;
    std::vector<double> mpi = MPI(res, tolerance2, b, x0, tau);
    std::vector<double> gaus = G_Z(res, tolerance2, b, x0);
    std::vector<double> gradient = fastest_gradient_descent(res, tolerance2, b, x0);

//std::vector<double> jac = Jacoby(res, tolerance2, b, x0);
    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(mpi[i], gradient[i], 1e-13);
    }
}

TEST(Grad, Gradient_second_test) {
    std::vector<element<double>> matrix_CSR = {
            {1, 5, 3},
            {4, 1, 17},
            {0, 25, 3},
            {9, 0, 7},
            {1, 1, 2},
            {1, 2, 28},
            {5, 1, 3},
            {6, 1, 189},
            {9, 2, 4}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res = CompressedMatrix(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0 = {1, 1, 1};

    double tau = 0.0001;
    double tolerance1 = 1e-12;
    double tolerance2 = 1e-15;
    std::vector<double> mpi = MPI(res, tolerance2, b, x0, tau);
    std::vector<double> gz = G_Z(res, tolerance2, b, x0);
    std::vector<double> Grad = fastest_gradient_descent(res, tolerance2, b, x0);
    //std::vector<double> jac = Jacoby(res, tolerance2, b, x0);
    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(mpi[i], Grad[i], 1e-13);
    }
}


