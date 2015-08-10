#include <zee.hpp>

#include "medium_grain/medium_grain.hpp"

#include <string>

using namespace Zee;

int main()
{
    std::string matrix = "steam3";

    ZeeLogInfo << "-- Starting IR example" << endLog;

    // Initialize the centralized base matrix from file
    auto baseMatrix = DSparseMatrix<double, int>("data/matrices/" + matrix  + ".mtx");
    baseMatrix.spy(matrix);
    auto& A = baseMatrix;

    vector<int> communicationVolumes;

    //-------------------------------------------------------------------------
    // MEDIUM GRAIN + IR
    MGPartitioner<decltype(baseMatrix)> mgPartitioner;
    mgPartitioner.partition(A);
    A.spy(matrix + "_initial_mg");

    communicationVolumes.push_back(A.communicationVolume());
    auto iter = 0;
    auto iters = 1000;
    while (!mgPartitioner.locallyOptimal() && iter < iters) {
        mgPartitioner.refine(A);
        communicationVolumes.push_back(A.communicationVolume());
        A.spy(matrix + "_mg_refine");
        ++iter;
    }
    ZeeLogVar(communicationVolumes);

    // the final matrix
    A.spy(matrix + "_mg_final");

    return 0;
}
