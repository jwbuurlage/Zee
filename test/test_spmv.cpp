#include <cstdint>
#include "zee.hpp"

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix with set fillin
    int32_t n = 1000;
    int32_t m = 1000;
    double fill_in = 0.1;

    DSparseMatrix<double> A(n, m);
     //A = rand(n, m, fill_in);

    // 4. Multiply A with some dense vector and store the result
//    DVector<double> v = rand(n);
//    DVector<double> u(n);
//
//    // compute spmv
//    SpMV<P_CPP>(A, v, u);

   return 0;
}
