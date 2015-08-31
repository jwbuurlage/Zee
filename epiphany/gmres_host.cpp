#include <string>
#include <cmath>

#include <unpain_cpp.hpp>
#include <zee.hpp>
#include "streaming_matrix.hpp"
#include "../examples/solvers/gmres.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = int;

    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    auto center = std::make_shared<Unpain::Center<TIdx>>(4);

    // We initialize the matrix with cyclic distribution
    std::string matrix = "steam3";
    auto A = DStreamingSparseMatrix<TVal, TIdx>(center, "../data/matrices/" + matrix + ".mtx", 4);

    auto b = DStreamingVector<TVal, TIdx>{center, A.getRows(), 1.0};
    b = A * b;

    // initial x is the zero vector
    auto x = DStreamingVector<TVal, TIdx>{center, A.getCols()};

    // Start GMRES
    GMRES::solve(A, // Matrix
            b,                  // RHS vector
            x,                  // resulting guess for x
            100,                // outer iterations
            100,                // inner iterations
            1e-6,               // tolerance level
            false);             // plot residuals

    DStreamingVector<TVal, TIdx> c{center, A.getRows()};
    DStreamingVector<TVal, TIdx> r{center, A.getRows()};
    c = A * x;
    r = b - c;
    ZeeLogVar(r.norm() / b.norm());

    return 0;
}
