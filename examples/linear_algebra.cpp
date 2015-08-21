#include <zee.hpp>

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = int;

    constexpr auto size = 10;

    // Vectors
    auto x = DVector<TVal, TIdx>{size, 1.0};
    auto y = DVector<TVal, TIdx>{size, 1.0};
    auto z = DVector<TVal, TIdx>{size, 1.0};

    z = x + y;
    ZeeLogVar(z);

    z = x + y + x;
    ZeeLogVar(z);

    z = (x + y) - x;
    ZeeLogVar(z);

    z = x + (y - x);
    ZeeLogVar(z);

    z = (x + y) + (x + y);
    ZeeLogVar(z);

    z = x * 4.5;
    ZeeLogVar(z);

    z = 4.5 * x;
    ZeeLogVar(z);

    z = x / 2.0;
    ZeeLogVar(z);

    // Matrices
    std::string matrix = "karate";
    auto A = DSparseMatrix<TVal, TIdx>("data/matrices/" + matrix + ".mtx", 4);
    auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
    auto u = DVector<TVal, TIdx>{A.getRows(), 1.0};

    u = A * v;
    ZeeLogVar(u);

    u = A * (v + v) + v;
    ZeeLogVar(u);

    // ...
}
