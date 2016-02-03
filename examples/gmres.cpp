#include <string>
#include <cmath>

#include <zee.hpp>
#include "solvers/gmres.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = unsigned int;

    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    // We initialize the matrix with cyclic distribution
    std::string matrix = "GD95_c";
    auto A = DSparseMatrix<TVal, TIdx>("data/matrices/" + matrix + ".mtx", 1);

    auto b = DVector<TVal, TIdx>{A.getRows(), 1.0};
    GreedyVectorPartitioner<decltype(A), decltype(b)> pVecs(A, b, b);
    pVecs.partition();
    pVecs.localizeMatrix();

    b = A * b;

    // initial x is the zero vector
    auto x = DVector<TVal, TIdx>{A.getCols()};

    // Start GMRES
    GMRES::solve<TVal, TIdx>(A, // Matrix
            b,                  // RHS vector
            x,                  // resulting guess for x
            1,                // outer iterations
            50,                // inner iterations
            1e-6,               // tolerance level
            true);             // plot residuals

    DVector<TVal, TIdx> c{A.getRows()};
    DVector<TVal, TIdx> r{A.getRows()};
    c = A * x;
    r = b - c;
    ZeeLogVar(r.norm() / b.norm());

    return 0;
}
