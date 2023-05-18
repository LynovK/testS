#file(GLOB_RECURSE all *.cpp *.h)
#add_library(slae STATIC ${all} "Compressed_spars_row/CSR.cpp" Tridiagonal/ThreeDiagonalSolver.h ../gtest/ThreeDiagTest.cpp)
#
#add_library(slae STATIC ${all} Tridiagonal/ThreeDiagonalSolver.h)

file(GLOB_RECURSE source *.h *.cpp)
add_library(slae INTERFACE ${source} Iteration_methods/Iteration_methods.h Compressed_spars_row/CompressedSorting.h Iteration_methods/ChebishevMPI.h Iteration_methods/Grad_spusk.h "rv_to_Kepler/rv_to_Kepler.h" Iteration_methods/Related_directions_method.h Iteration_methods/Symmetric_G_Z.h Vector.h Vector.h Vector.h operator/Operators.h Tridiagonal/five_diag.h)
