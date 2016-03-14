#include <zee.hpp>

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

    JWAssert(A.getRows() == A.getCols());

    auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
    auto u = DVector<TVal, TIdx>{A.getRows()};

    MGPartitioner<decltype(A)> partitioner;
    partitioner.partition(A);

    JWLogVar(A.communicationVolume());
    JWLogVar(A.loadImbalance());

    GreedyVectorPartitioner<decltype(A), decltype(v)> pVecs(A, v, u);
    pVecs.partition();

    JWLogVar(v.getOwners());
    JWLogVar(u.getOwners());

    pVecs.localizeMatrix();

    JWLogVar(u);
    u = A * v;
    JWLogVar(u);

    return 0;
}
