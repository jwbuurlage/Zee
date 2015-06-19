#include <zee.hpp>

#include "pulp/pulp.hpp"
#include "medium_grain/medium_grain.hpp"

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

    DSparseMatrix<double, int> A = fromMM<double, int>("data/matrices/karate.mtx", 1);

    A.spy();

    logError("phony");

    //PulpPartitioner<decltype(A)> partitioner;
    MGPartitioner<decltype(A)> partitioner;
    partitioner.initialize(A);

    auto& B = partitioner.partition(A);
    B.spy();

    Zee::CyclicPartitioner<decltype(A)> cyclicPartitioner(2, CyclicType::column);
    auto& C = cyclicPartitioner.partition(A);
    C.spy();

    return 0;
}
