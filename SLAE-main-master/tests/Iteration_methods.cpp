#include <gtest/gtest.h>
#include <vector>
#include "../src/Compressed_spars_row/CompressedMatrix.h"
#include "../src/Compressed_spars_row/CompressedSorting.h"
#include "../src/Iteration_methods/Iteration_methods.h"
#include "../src/Iteration_methods/Symmetric_G_Z.h"
#include <fstream>


TEST(task_4, mpi_test) {
    std::vector<element<double>> matrix_CSR = {
            {6, 4, 11},
            {2, 9, 13},
            {0, 2, 3},
            {13, 1, 1},
            {10, 44, 2},
            {5, 8, 28},
            {2, 0, 3},
            {1, 2, 45},
            {2, 2, 69}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 2, 3};
    std::vector<double> x0 = {1, 1, 1};
    double tau = 0.0001;
    double tolerance1 = 1e-12;
    double tolerance2 = 1e-15;
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau);

    for (auto i: mpi) {
        std::cout << i << std::endl;
    }
}

TEST(task_4, SOR_test) {
    std::vector<element<double>> matrix_CSR = {
            {6, 3, 10},
            {5, 1, 2},
            {2, 2, 1},
            {1, 1, 2},
            {2, 4, 10},
            {4, 2, 3},
            {2, 5, 1},
            {7, 1, 3},
            {2, 3, 10}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res(matrix_CSR, 3, 3);

    std::vector<double> b = {1, 1, 1};
    std::vector<double> x0 = {0, 0, 0};
    double tau = 0.0001;
    double tolerance1 = 1e-12;
    double omega = 0.5;
    std::vector<double> sor = SOR(res, tolerance1, omega, b, x0);
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau);
    for (auto i: sor) {
        std::cout << i << std::endl;
    }
    for (auto i = 0; i < mpi.size(); ++i) {
        ASSERT_NEAR(sor[i], mpi[i], 1e-12);
    }
}

TEST(task_4, Gaus_Zeidel_firsttest) {
    std::vector<element<double>> matrix_CSR = {
            {4, 1, 12},
            {1, 4, 17},
            {8, 5, 3},
            {14, 0, 17},
            {5, 1, 13},
            {5, 2, 28},
            {2, 0, 3},
            {4, 1, 1},
            {2, 6, 9}
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
            {3, 1, 10},
            {3, 6, 2},
            {0, 0, 1},
            {8, 0, 2},
            {4, 1, 10},
            {3, 24, 3},
            {2, 7, 1},
            {2, 6, 3},
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
    std::vector<double> gaus = G_Z(res, tolerance2, b, x0);

    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(mpi[i], gaus[i], 1e-10);
    }
}


TEST(sym_G_Z, symGZ_first) {
    std::vector<element<double>> matrix_CSR = {
            {6, 5, 10},
            {2, 1, 2},
            {3, 2, 1},
            {1, 7, 2},
            {3, 3, 1},
            {4, 4, 3},
            {8, 0, 1},
            {2, 2, 2},
            {2, 2, 4}
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

    for (long i = 0; i < 3; i++) {
        ASSERT_NEAR(gaus[i], jacoby[i], 1e-10);
    }
}


TEST(difuri, Difur_dz_test) {
    std::vector<element<double>> matrix_CSR = {
            {0, 0, 1},
            {0, 1, 0},
            {0, 2, 1},
            {0, 3, 1},
            {1, 0, 0},
            {1, 1, 1},
            {1, 2, 2},
            {1, 3, -2},
            {2, 0, 1},
            {2, 1, 1},
            {2, 2, M_E * M_E},
            {2, 3, 1 / (M_E * M_E)},
            {3, 0, 0},
            {3, 1, 1},
            {3, 2, 2 * M_E * M_E},
            {3, 3, -2 / (M_E * M_E)}
    };
    sort_me_plz(matrix_CSR);
    CompressedMatrix<double> res(matrix_CSR, 4, 4);

    std::vector<double> b = {0, 0, 0.25 * (M_E * M_E - 3), 0.5 * (M_E * M_E - 1)};
    std::vector<double> x0 = {0, 0, 0, 0};
    double tau = 0.01;
    double tolerance1 = 1e-12;
    double omega = 0.5;
    std::vector<double> sor = SOR(res, tolerance1, omega, b, x0);
    std::vector<double> mpi = MPI(res, tolerance1, b, x0, tau);
    for (auto i: sor) {
        std::cout << i << std::endl;
    }
    for (auto i = 0; i < mpi.size(); ++i) {
        ASSERT_NEAR(sor[i], mpi[i], 1e-12);
    }
    for (auto i: mpi) {
        std::cout << i << std::endl;
    }
}
