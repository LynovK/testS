#include <iostream>
#include <gtest/gtest.h>
#include "../src/Iteration_methods/ChebishevMPI.h"

TEST(task_4, task_test_twotwo) {
    std::vector<element<double>> matrix_CSR = {
            {0, 1, 2},
            {1, 2, 3},
            {2, 3, 4},
            {3, 4, 5},
            {4, 5, 6},
            {5, 6, 7},
            {6, 7, 8},
            {7, 8, 9},
            {8, 9, 10},
            {9, 10, 11}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res(matrix_CSR, 10, 10);
    std::vector<double> tau = Chebishev_solver<double>(16, 4, 1, 10);

    std::vector<double> solve = {5,
                                 1,
                                 0.667,
                                 1.75,
                                 1.6,
                                 7.16667,
                                 0.2857,
                                 0.125,
                                 0.556,
                                 0.6};

    std::vector<double> b = {5, 2, 2, 7, 8, 43, 2, 1, 5, 6};
    std::vector<double> x0 = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    double tau1 = 0.001;
    double tolerance1 = 1e-12;
    double tolerance2 = 1e-15;
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau1);
    std::vector<double> mpi_cheb = MPI_chebishev(res, tolerance1, b, x0, tau);


    for(int i = 0; i < solve.size(); ++i){
        ASSERT_NEAR(mpi_cheb[i], mpi[i], 1e-12);
    }
}



    double lambda_max = 24.0;
    double tau1 = 1.8/lambda_max;

    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res(matrix_CSR, 4, 4);
    std::vector<double> tau = Chebishev_solver<double>(8, 3, 16, 24);

    std::vector<double> b = {6, 6, 6, 6};
    std::vector<double> b_pop = {3, 3, 3, 3};
    std::vector<double> x0 = {1, 1, 1, 1};

    double tolerance1 = 1e-13;
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau1);
    std::vector<double> mpi_cheb = MPI_chebishev(res, tolerance1, b, x0, tau);

    for(int i = 0; i < mpi.size(); ++i){
        ASSERT_NEAR(mpi_cheb[i], mpi[i], 1e-13);
    }
}
