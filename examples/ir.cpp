#include <zee.hpp>

#include "medium_grain/medium_grain.hpp"

#include <string>

using namespace Zee;

int main()
{
    std::string matrix = "karate";

    ZeeLogInfo << "-- Starting IR example" << endLog;

    // Initialize the centralized base matrix from file
    DSparseMatrix<double, int> baseMatrix =
        fromMatrixMarket<double, int>("data/matrices/" + matrix  + ".mtx", 1);
    baseMatrix.spy("karate");
    auto& A = baseMatrix;

    vector<int> communicationVolumes;

    //-------------------------------------------------------------------------
    // MEDIUM GRAIN + IR
    MGPartitioner<decltype(baseMatrix)> mgPartitioner;
    mgPartitioner.partition(A);
    A.spy(matrix + "_initial_mg");

    communicationVolumes.push_back(A.communicationVolume());
    while (!mgPartitioner.locallyOptimal()) {
        mgPartitioner.refine(A);
        communicationVolumes.push_back(A.communicationVolume());
        A.spy(matrix + "_initial_mg_refine");
    }
    ZeeLogVar(communicationVolumes);

    // the final matrix
    A.spy(matrix + "_initial_mg_final");

    return 0;
}
