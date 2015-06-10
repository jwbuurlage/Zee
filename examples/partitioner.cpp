#include <zee.hpp>

#include "pulp/pulp.hpp"

#include <cstdint>
#include <iostream>
using std::cout;
using std::endl;

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix with set fillin
    uint32_t n = 20;
    uint32_t m = 20;
    uint32_t p = 3;
    double fill_in = 0.4;
    DSparseMatrix<double> A = rand(n, m, p, fill_in);

    cout << A.loadImbalance() << endl;
    cout << A.communicationVolume() << endl;

    // 2a. Initialize a cyclic partitioner
    // pulp is 'iterative', used to refine
    A.spy();
    PulpPartitioner<DSparseMatrix<double>> pulpPart;
    for (int i = 0; i < 1000; ++i) {
        pulpPart.refine(A);
    }
    A.spy();

    return 0;
}
