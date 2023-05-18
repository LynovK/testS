//
// Created by perseverance on 09.04.2023.
//

#include <gtest/gtest.h>
#include <vector>
#include "../src/Compressed_spars_row/CompressedMatrix.h"
#include "../src/Compressed_spars_row/CompressedSorting.h"
#include "../src/Iteration_methods/Iteration_methods.h"
#include "../src/Iteration_methods/ChebishevMPI.h"
#include "../src/Tridiagonal/five_diag.h"

TEST(KR2, first_task_test) {
    constexpr unsigned int n = 289;
    constexpr double b = 9;
    constexpr double a = 4;
    constexpr unsigned int L = 17;

    five_diag<double> matrix_(n, a, b, L);

    sort_me_plz(matrix_.matrix);
    CompressedMatrix<double> res(matrix_.matrix, n, n);
    std::vector<double> d(n, 4);
    std::vector<double> x0(n, 0);
    double tolerance1 = 1e-13;

    // тау пункт 1
    double k = 2 * (b + 2 * a * cos(M_PI / (L + 1)));
    double tau = 1 / k;
    std::vector<double> grad = MPI(res, tolerance1, d, x0, tau);
    //

    // град спуск с оптимальным тау (пункт 2)
    double tau_opt = 2 / (matrix_.lambda_max() + matrix_.lambda_min());
    std::vector<double> grad_opt = MPI(res, tolerance1, d, x0, tau_opt);
    //

    // ускорение чебышева (пункт 3)
    std::vector<double> tauu = Chebishev_solver<double>(8, 3, matrix_.lambda_min(), matrix_.lambda_max());
    std::vector<double> mpi_cheb = MPI_chebishev(res, tolerance1, d, x0, tauu);

    for (auto i = 0; i < grad.size(); ++i) {
        ASSERT_NEAR(grad[i], grad_opt[i], 1e-13);
    }
}

TEST(KR2, second_task_test) {
    std::vector<element<double>> matrix_CSR = {
            {0, 0, 16},
            {1, 1, 18.0},
            {2, 2, 21.0},
            {3, 3, 24.0},
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res(matrix_CSR, 4, 4);
    std::vector<double> b = {6, 6, 6, 6};
    std::vector<double> x0 = {0, 0, 0, 0};
    double tolerance1 = 1e-13;

//град спуск пункт 1
    double lambda_max = 24.0;
    double lambda_min = 16;

    double tau1 = 1.8 / lambda_max;

    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau1);

//град спуск пункт 2 оптимальное тау
    double lambda_opt = 2 / (lambda_max + lambda_min);
//std::vector<double> mpi_opt = MPI(res, tolerance1, b, x0, lambda_opt);

//гаиск град спуск пункт 3
//std::vector<double> naisk_grad = fastest_gradient_descent(res, tolerance1, b, x0);

//град спуск чебышев пункт 4
    std::vector<double> tau = Chebishev_solver<double>(8, 3, lambda_min, lambda_max);
    std::vector<double> mpi_cheb = MPI_chebishev(res, tolerance1, b, x0, tau);

//сопр градиенты
//std::vector<double> CG = Related_directions(res, tolerance1, b, x0);

/*    for(int i = 0; i < mpi.size(); ++i){
        ASSERT_NEAR(mpi_cheb[i], mpi[i], 1e-13);
    }*/
}