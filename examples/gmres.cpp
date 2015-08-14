#include <string>
#include <cmath>

#include <unpain_cpp.hpp>
#include <zee.hpp>
#include "solvers/gmres.hpp"

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = int;

    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    auto center = std::make_shared<Unpain::Center<TIdx>>(4);

    // We initialize the matrix with cyclic distribution
    std::string matrix = "karate";
    auto A = DSparseMatrix<TVal, TIdx>(center, "data/matrices/" + matrix + ".mtx", 4);

    auto b = DVector<TVal, TIdx>{center, A.getRows(), 1.0};
    b = A * b;

    // initial x is the zero vector
    auto x = DVector<TVal, TIdx>{center, A.getCols()};

    // Start GMRES
    GMRES::solve<TVal, TIdx>(A, // Matrix
            b,                  // RHS vector
            x,                  // resulting guess for x
            100,                // outer iterations
            100,                // inner iterations
            1e-6);              // tolerance level

    DVector<TVal, TIdx> c{center, A.getRows()};
    DVector<TVal, TIdx> r{center, A.getRows()};
    c = A * x;
    r = b - c;
    ZeeLogVar(r.norm() / b.norm());

    return 0;
}
