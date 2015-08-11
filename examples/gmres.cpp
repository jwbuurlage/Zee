#include <string>

#include <unpain_cpp.hpp>
#include <zee.hpp>

#include "medium_grain/medium_grain.hpp"

using namespace Zee;

namespace GMRES
{

template<typename TVal, typename TIdx>
void solve(const DSparseMatrix<TVal, TIdx>& A,
        const DVector<TVal, TIdx>& b,
        DVector<TVal, TIdx>& x)
{
    auto& center = A.getCenter();

    auto r = b;
    DVector<TVal, TIdx> e{center, A.rows()};

    auto maxIterations = 10;

    for (auto iteration = 0; iteration < maxIterations; ++iteration)
    {

        ZeeLogInfo << "|| b - A x || = : " << r.norm() << endLog;
    }
}

} // namespace GMRES

int main()
{
    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    auto b = Benchmark("GMRES");

    auto center = std::make_shared<Unpain::Center<int>>(2);
    std::string matrix = "karate";

    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx");
    auto b = DVector<double, int>{center, A.getCols(), 1.0};

    MGPartitioner<decltype(A)> partitioner;
    partitioner.partition(A);
    for (auto iter = 0; iter < 100; ++iter) {
        partitioner.refine(A);
    }
    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());

    // initial x is the zero vector
    auto x = DVector<double, int>{center, A.getRows()};

    // Start GMRES
    GMRES::solve<double, int>(A, b, x);

//    auto r = b - A * x;
//    ZeeLogInfo << "|| b - A x || = : " << r.norm() << endLog;

    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());


    return 0;
}
