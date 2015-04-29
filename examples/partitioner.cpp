#include <zee.h>
#include <pulp/pulp.hpp>

#include <cstdint>
#include <iostream>
using std::cout;
using std::endl;

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix with set fillin
    uint32_t n = 40;
    uint32_t m = 40;
    uint32_t p = 4;
    double density = 0.1;
    DSparseMatrix<double> A = rand(n, m, density, p);

    cout << A.loadImbalance() << endl;
    cout << A.communicationVolume() << endl;

    // 2a. Initialize a cyclic partitioner
    // pulp is 'iterative', used to refine
    A.spy();
    PulpPartioner<double> pulpPart;
    pulpPart.refine(&A);
    A.spy();

    // 2b. bisectioning partitioner (n -> 2n)
    DSparseMatrix<double> B = rand(n, m, fill_in, p);

    B.spy();
    BisectionPartitioner<double> bisectPart;
    bisectPart.refine(&B);
    B.spy();

    // 4. Multiply A with some dense vector and store the result
    DDenseMatrix<double> v = rand(n, 1);
    DDenseMatrix<double> u;
    u = A * v;

    return 0;
}
