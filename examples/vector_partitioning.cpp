#include <zee.hpp>

#include "partitioners/pulp.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = uint32_t;
    std::string matrix = "karate";

    // Initialize the centralized base matrix from file
    DSparseMatrix<TVal, TIdx> baseMatrix{"data/matrices/" + matrix  + ".mtx", 1};
    auto& A = baseMatrix;
    auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
    auto u = DVector<TVal, TIdx>{A.getCols()};

    PulpPartitioner<decltype(A)> pA(A);
    pA.refineWithIterations(1000);

    GreedyVectorPartitioner<decltype(v)> pVecs(A, v, u);
    pVecs.partition();

    A.localizeIndices(v, u);

    return 0;
}
