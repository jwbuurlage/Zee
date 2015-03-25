#include <cstdint>
#include "zee.hpp"

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix with set fillin
    uint32_t n = 1000;
    uint32_t m = 1000;
    uint32_t p = 4;
    double fill_in = 0.1;

    DSparseMatrix<double> A = rand(n, m, p, fill_in);

    // A.debugOutput();

    // for now vectors are not distributed
    DVector<double> v = rand(n);
    DVector<double> u(n);

    // 4. Multiply A with some dense vector and store the result
    //SpMV<P_CPP>(A, v, u);

   return 0;
}