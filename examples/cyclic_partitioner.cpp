#include <zee.hpp>

#include <cstdint>
#include <iostream>
using std::cout;
using std::endl;

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix on a single processor
    uint32_t n = 20;
    uint32_t m = 20;
    uint32_t p = 1;
    double fill_in = 0.4;
    DSparseMatrix<double> A = rand(n, m, p, fill_in);

    // Spy it
    A.spy();

    Factory<CyclicPartitioner<decltype(A)>> f;

    {
        auto cp = f.make();
        cp->setProcs(4);
        cp->partition(A);
    }

    // Spy it again
    A.spy();

    return 0;
}
