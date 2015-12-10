#include <zee.hpp>
#include "partitioners/pulp.hpp"

using namespace Zee;

int main()
{
    using TVal = float;
    using TIdx = uint32_t;
    std::string matrix = "mesh1e1";
    TIdx procs = 1;

    // Initialize the matrix from file
    DSparseMatrix<TVal, TIdx> A{
        "data/matrices/" + matrix  + ".mtx",
        procs
    };

    ZeeAssert(A.getRows() == A.getCols());

    auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
    auto u = DVector<TVal, TIdx>{A.getRows()};

    PulpPartitioner<decltype(A)> pA(A);
    pA.refineWithIterations(1000);

    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());

    GreedyVectorPartitioner<decltype(A), decltype(v)> pVecs(A, v, u);
    pVecs.partition();

    ZeeLogVar(v.getOwners());
    ZeeLogVar(u.getOwners());

    pVecs.localizeMatrix();

    ZeeLogVar(u);
    u = A * v;
    ZeeLogVar(u);

    return 0;
}
