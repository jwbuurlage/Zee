#include <string>
#include <cmath>

#include <zee.hpp>
#include "streaming_matrix.hpp"
#include "../../examples/solvers/gmres.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = int;

    ZeeLogInfo << "-- Starting Epiphany example" << endLog;

    // We initialize the matrix with cyclic distribution
    std::string matrix = "steam3";
    auto A = DStreamingSparseMatrix<TVal, TIdx>("../../data/matrices/" + matrix + ".mtx", 4);
    DStreamingVector<TVal, TIdx> x(A.getCols(), 1.0);
    DStreamingVector<TVal, TIdx> y(A.getRows(), 1.0);
    y = A * x;

    return 0;
}
