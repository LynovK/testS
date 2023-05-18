file(GLOB_RECURSE source *.h *.cpp)
add_library(slae INTERFACE ${source} Iteration_methods/Iteration_methods.h Compressed_spars_row/CompressedSorting.h Iteration_methods/ChebishevMPI.h Iteration_methods/Grad_spusk.h  Iteration_methods/Related_directions_method.h Iteration_methods/Symmetric_G_Z.h Vector.h Vector.h Vector.h operator/Operators.h Tridiagonal/five_diag.h)
