#include <gtest/gtest.h>
#include <vector>
#include "../src/Compressed_spars_row/CompressedMatrix.h"
#include "../src/Compressed_spars_row/CompressedSorting.h"
#include "../src/Iteration_methods/Iteration_methods.h"
#include "../src/Iteration_methods/Symmetric_G_Z.h"
#include "../src/Iteration_methods/ChebishevMPI.h"
#include "../src/Iteration_methods/Grad_spusk.h"
#include <fstream>

TEST(task_4, Gaus_Zeidel_firsttest) {
    std::vector<element<double>> matrix_CSR = {
            {0, 0, 12},
            {0, 1, 17},
            {0, 2, 3},
            {1, 0, 17},
            {1, 1, 15825},
            {1, 2, 28},
            {2, 0, 3},
            {2, 1, 28},
            {2, 2, 238}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res = CompressedMatrix(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0 = {1, 1, 1};

    double tau = 0.0001;
    double tolerance1 = 1e-12;
    double tolerance2 = 1e-15;
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau);
    std::vector<double> gz = G_Z(res, tolerance2, b, x0);
    //std::vector<double> jac = Jacoby(res, tolerance2, b, x0);
    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(mpi[i], gz[i], 1e-13);
    }
}

TEST(task_4, Gaus_Zeidel_secondtest) {
    std::vector<element<double>> matrix_CSR = {
            {0, 0, 10},
            {0, 1, 2},
            {0, 2, 1},
            {1, 0, 2},
            {1, 1, 10},
            {1, 2, 3},
            {2, 0, 1},
            {2, 1, 3},
            {2, 2, 10}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res = CompressedMatrix(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 1, 1};
    std::vector<double> x0 = {0, 0, 0};

    double tau = 0.0001;
    double tolerance1 = 1e-10;
    double tolerance2 = 1e-10;
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau);
    //std::vector<double> gz = G_Z(res, tolerance2, b, x0);
    std::vector<double> gaus = G_Z(res, tolerance2, b, x0);

    //std::vector<double> jac = Jacoby(res, tolerance2, b, x0);
    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(mpi[i], gaus[i], 1e-10);
    }
}

//symmetric tests

TEST(sym_G_Z, symGZfirst) {
    std::vector<element<double>> matrix_CSR = {
            {0, 0, 10},
            {0, 1, 2},
            {0, 2, 1},
            {1, 0, 2},
            {1, 1, 10},
            {1, 2, 3},
            {2, 0, 1},
            {2, 1, 3},
            {2, 2, 10}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res = CompressedMatrix(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 1, 1};
    std::vector<double> x0 = {0, 0, 0};

    double tau = 0.001;
    double tolerance1 = 1e-10;
    double tolerance2 = 1e-10;
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau);
    //std::vector<double> gz = G_Z(res, tolerance2, b, x0);
    std::vector<double> gaus = G_Z(res, tolerance2, b, x0);
    std::vector<double> sym_gaus = Symmentric_G_Z(res, tolerance2, b, x0);
    std::vector<double> jacoby = Jacoby(res, tolerance2, b, x0);

    //std::vector<double> jac = Jacoby(res, tolerance2, b, x0);
/*    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(mpi[i], gaus[i], 1e-10);
    }*/
    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(gaus[i], jacoby[i], 1e-10);
    }
}
