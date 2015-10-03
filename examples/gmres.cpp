#include <string>
#include <cmath>

#include <zee.hpp>
#include "solvers/gmres.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = int;

    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    // We initialize the matrix with cyclic distribution
    std::string matrix = "karate";
    auto A = DSparseMatrix<TVal, TIdx>("data/matrices/" + matrix + ".mtx", 4);

    auto b = DVector<TVal, TIdx>{A.getRows(), 1.0};
    b = A * b;

    // initial x is the zero vector
    auto x = DVector<TVal, TIdx>{A.getCols()};

    // Start GMRES
    GMRES::solve<TVal, TIdx>(A, // Matrix
            b,                  // RHS vector
            x,                  // resulting guess for x
            100,                // outer iterations
            100,                // inner iterations
            1e-6,               // tolerance level
            false);             // plot residuals

    DVector<TVal, TIdx> c{A.getRows()};
    DVector<TVal, TIdx> r{A.getRows()};
    c = A * x;
    r = b - c;
    ZeeLogVar(r.norm() / b.norm());

    return 0;
}
