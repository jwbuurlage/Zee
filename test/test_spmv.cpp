#include <cstdint>

#include <iostream>

#include "zee.hpp"

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix with set fillin
    uint32_t n = 30;
    uint32_t m = 30;
    uint32_t p = 4;
    double fill_in = 0.1;

    DSparseMatrix<double> A = rand(n, m, p, fill_in);
    // auto A = eye(n, p);
    A.spy();

    // for now vectors are not distributed
    DVector<double> v = rand(n);
    DVector<double> u = zeros(n);

    // 4. Multiply A with some dense vector and store the result
    spmv_cpp<double, uint32_t>(A, v, u);

    return 0;
}
