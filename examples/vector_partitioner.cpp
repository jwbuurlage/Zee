#include <zee.hpp>

#include "partitioners/greedy_vector.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = uint32_t;
    std::string matrix = "karate";
    TIdx procs = 2;

    // Initialize the matrix from file
    DSparseMatrix<TVal, TIdx> A{
        "data/matrices/" + matrix  + ".mtx",
        procs
    };


    auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
    auto u = DVector<TVal, TIdx>{A.getCols()};

//    PulpPartitioner<decltype(A)> pA(A);
//    pA.refineWithIterations(1000);
    MGPartitioner<decltype(A)> partitioner;
    partitioner.partition(A);

    GreedyVectorPartitioner<decltype(A), decltype(v)> pVecs(A, v, u);
    pVecs.partition();

    ZeeLogVar(v.getOwners());
    ZeeLogVar(u.getOwners());

    pVecs.localizeMatrix();

    return 0;
}
