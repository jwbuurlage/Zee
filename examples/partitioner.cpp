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
    for (int i = 0; i < 1000; ++i) {
        pulpPart.refine(A);
    }
    A.spy();

    return 0;
}
