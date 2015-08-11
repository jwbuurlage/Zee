#include <string>

#include <unpain_cpp.hpp>
#include <zee.hpp>

#include "medium_grain/medium_grain.hpp"

using namespace Zee;

int main()
{
    ZeeLogInfo << "-- Starting SpMV example" << endLog;

    auto b = Benchmark("SpMV");
    b.phase("load matrix and vectors");

    auto center = std::make_shared<Unpain::Center<int>>(2);
    std::string matrix = "fpga_dcop_05";

    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx");
    auto v = DVector<double, int>{center, A.getCols(), 1.0};
    auto u = DVector<double, int>{center, A.getRows()};

    b.phase("partitioning");

    MGPartitioner<decltype(A)> partitioner;
    partitioner.partition(A);
    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());

    b.phase("spmv (1)");
    for (int i = 0; i < 10; ++i) {
        u = A * v;
    }

    b.phase("refinement");
    for (auto iter = 0; iter < 100; ++iter) {
        partitioner.refine(A);
    }
    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());
    A.spy("fpga_refined");

    b.phase("spmv (2)");
    for (int i = 0; i < 10; ++i) {
        u = A * v;
    }

    return 0;
}
