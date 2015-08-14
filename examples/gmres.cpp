#include <string>
#include <cmath>

#include <unpain_cpp.hpp>
#include <zee.hpp>
#include "solvers/gmres.hpp"

using namespace Zee;

int main()
{
    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    auto center = std::make_shared<Unpain::Center<int>>(4);

    // We initialize the matrix with cyclic distribution
    std::string matrix = "karate";
    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx", 4);

    auto b = DVector<double, int>{center, A.getRows(), 1.0};
    b = A * b;

    // initial x is the zero vector
    auto x = DVector<double, int>{center, A.getCols()};

    // Start GMRES
    GMRES::solve<double, int>(A, // Matrix
            b,                   // RHS vector
            x,                   // resulting guess for x
            100,                  // outer iterations
            100,                  // inner iterations
            1e-6);               // tolerance level

    DVector<double, int> c{center, A.getRows()};
    DVector<double, int> r{center, A.getRows()};
    c = A * x;
    r = b - c;
    ZeeLogVar(r.norm() / b.norm());

    return 0;
}
