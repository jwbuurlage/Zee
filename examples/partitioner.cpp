#include <zee.hpp>

#include "pulp/pulp.hpp"
#include "medium_grain/medium_grain.hpp"
#include "multi_level/multi_level.hpp"

#include <string>
#include <cstdint>
#include <iostream>
using std::cout;
using std::endl;

using namespace Zee;

int main()
{
    std::string matrix = "karate";

    ZeeLogInfo << "-- Starting partitioning example" << endLog;

    // Initialize the centralized base matrix from file
    DSparseMatrix<double, int> baseMatrix =
        fromMatrixMarket<double, int>("data/matrices/" + matrix  + ".mtx", 1);
    baseMatrix.spy("karate");
    auto& A = baseMatrix;

    //-------------------------------------------------------------------------
    // MEDIUM GRAIN
    MGPartitioner<decltype(baseMatrix)> mgPartitioner;
    mgPartitioner.partition(A);
    A.spy(matrix + "_mg");

    ZeeLogInfo << "MG: \t" << A.communicationVolume() << endLog;
    while (!mgPartitioner.locallyOptimal()) {
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
