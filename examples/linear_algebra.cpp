#include <zee.hpp>

using namespace Zee;

int main()
{
    using TVal = double;
    using TIdx = unsigned int;

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

    // Sparse Matrices
    std::string matrix = "karate";
    DSparseMatrix<TVal, TIdx> A("data/matrices/" + matrix + ".mtx", 4);
    DVector<TVal, TIdx> v{A.getCols(), 1.0};
    DVector<TVal, TIdx> u{A.getRows(), 1.0};

    u = A * v;
    ZeeLogVar(u);

    u = A * (v + v) + v;
    ZeeLogVar(u);

    // Dense matrices
    DMatrix<TVal, TIdx> D("data/matrices/dense_example.mtx");
    DMatrix<TVal, TIdx> E("data/matrices/dense_example.mtx");
    E.transpose();
    DMatrix<TVal, TIdx> F(D.getRows(), D.getRows());

    F = D * E;
    ZeeLogVar(D);
    ZeeLogVar(E);
    ZeeLogVar(F);
}
