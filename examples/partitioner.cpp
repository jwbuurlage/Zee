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

    ZeeInfoLog << "Starting partitioning example" << endLog;

    DSparseMatrix<double, int> A = fromMatrixMarket<double, int>("data/matrices/" + matrix  + ".mtx", 1);

    A.spy("karate");

    MGPartitioner<decltype(A)> partitioner;
    partitioner.initialize(A);
    auto& B = partitioner.partition(A);
    B.spy(matrix + "_mg");

    ZeeInfoLog << "MG: \t" << B.communicationVolume() << endLog;

    Zee::CyclicPartitioner<decltype(A)> cyclicPartitioner(8, CyclicType::column);
    auto& C = cyclicPartitioner.partition(A);
    C.spy(matrix + "_cyclic");

    ZeeInfoLog << "Cyclic: \t" << C.communicationVolume() << endLog;

    MultiLevelOneD<decltype(A)> mlPart{};
    mlPart.initialize(A);
    auto& D = mlPart.partition(A);
    C.spy(matrix + "_ml");

    ZeeInfoLog << "ML: \t" << D.communicationVolume() << ", " << D.loadImbalance() << endLog;

    return 0;
}
