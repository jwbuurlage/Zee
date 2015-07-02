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
    // 1. Initialize a random n*m sparse matrix with set fillin
    //uint32_t n = 100;
    //uint32_t m = 100;
    //uint32_t p = 4;
    //double fill_in = 0.4;
    //DSparseMatrix<double> A = rand(n, m, p, fill_in);
    
    std::string matrix = "karate";

    ZeeInfoLog << "Starting partitioner" << endLog;

    DSparseMatrix<double, int> A = fromMatrixMarket<double, int>("data/matrices/" + matrix  + ".mtx", 1);

    A.spy("karate");

    //PulpPartitioner<decltype(A)> partitioner;
    MGPartitioner<decltype(A)> partitioner;
    partitioner.initialize(A);
    auto& B = partitioner.partition(A);
    B.spy(matrix + "_mg");

    Zee::CyclicPartitioner<decltype(A)> cyclicPartitioner(8, CyclicType::column);
    auto& C = cyclicPartitioner.partition(A);
    C.spy(matrix + "_cyclic");

    ZeeInfoLog << "Cyclic: \t" << C.communicationVolume() << endLog;

    MultiLevelOneD<decltype(A)> mlPart{};
    mlPart.initialize(A);
    auto& D = mlPart.partition(A);
    C.spy(matrix + "_ml");

    ZeeInfoLog << "ML: \t" << D.communicationVolume() << endLog;


    return 0;
}
