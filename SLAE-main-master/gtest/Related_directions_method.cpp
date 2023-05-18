#include <gtest/gtest.h>
#include "../src/Iteration_methods/Related_directions_method.h"
#include "../src/Iteration_methods/ChebishevMPI.h"

TEST(Grad_napr, grad_napr_first){
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

    double tau_u = 0.001;
    double tolerance1 = 1e-12;
    double tolerance2 = 1e-15;
    std::vector<double> mpi_cheb = MPI_chebishev(res, tolerance2, b, x0, tau);
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau_u);
    std::vector<double> Grad_napr = Related_directions(res, tolerance2, b, x0);


    for(auto i = 0; i < solve.size(); ++i){
        ASSERT_NEAR(mpi_cheb[i], Grad_napr[i], 1e-15);
    } //макс точность 1e-15
}
