#include <zee.hpp>

#include "partitioners/pulp.hpp"
#include "partitioners/medium_grain.hpp"
#include "partitioners/multi_level.hpp"

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

    std::string matrix = "karate";

    ZeeLogInfo << "-- Starting partitioning example" << endLog;

    // Initialize the centralized base matrix from file
    DSparseMatrix<TVal, TIdx> baseMatrix{"data/matrices/" + matrix  + ".mtx", 1};
    baseMatrix.spy("karate");
    auto& A = baseMatrix;

    //-------------------------------------------------------------------------
    // MEDIUM GRAIN
    MGPartitioner<decltype(baseMatrix)> mgPartitioner;
    mgPartitioner.partition(A);
    A.spy(matrix + "_mg");

    ZeeLogInfo << "MG: \t" << A.communicationVolume() << endLog;
    while (!mgPartitioner.locallyOptimal()) {
        ZeeLogVar(A.communicationVolume());
        mgPartitioner.refine(A);
    }

    //-------------------------------------------------------------------------
    // 1D MULTILEVEL COLUMN
    MultiLevelOneD<decltype(baseMatrix)> mlPart{};
    mlPart.initialize(A);
    auto& D = mlPart.partition(A);
    D.spy(matrix + "_ml");

    ZeeLogInfo << "ML: \t" << D.communicationVolume() << ", " <<
        D.loadImbalance() << endLog;

    ZeeLogInfo << "-- Ending partitioning example" << endLog;

    return 0;
}
