#include <zee.hpp>

#include "partitioners/pulp.hpp"

#include <string>
#include <cstdint>
#include <iostream>
using std::cout;
using std::endl;

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = uint32_t;

    std::string matrix = "sparse_example";

    // initialize the centralized base matrix from file
    DSparseMatrix<TVal, TIdx> baseMatrix{"data/matrices/" + matrix  + ".mtx", 1};
    baseMatrix.spy(matrix);
    auto& A = baseMatrix;

    // medium grain
    MGPartitioner<decltype(baseMatrix)> mgPartitioner;
    mgPartitioner.partition(A);
    A.spy(matrix + "_mg");

    JWLogInfo << "MG: \t" << A.communicationVolume() << endLog;
    while (!mgPartitioner.locallyOptimal()) {
        JWLogVar(A.communicationVolume());
        mgPartitioner.refine(A);
    }

    // 1d multi-level by column
    MultiLevelOneD<decltype(baseMatrix)> mlPart{};
    mlPart.initialize(A);
    auto& D = mlPart.partition(A);
    D.spy(matrix + "_ml");

    JWLogInfo << "ML: \t" << D.communicationVolume() << ", " <<
        D.loadImbalance() << endLog;

    return 0;
}
