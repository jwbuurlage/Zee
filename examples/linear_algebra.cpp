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
    JWLogVar(z);

    z = x + y + x;
    JWLogVar(z);

    z = (x + y) - x;
    JWLogVar(z);

    z = x + (y - x);
    JWLogVar(z);

    z = (x + y) + (x + y);
    JWLogVar(z);

    z = x * 4.5;
    JWLogVar(z);

    z = 4.5 * x;
    JWLogVar(z);

    z = x / 2.0;
    JWLogVar(z);

    // Sparse Matrices
    DSparseMatrix<TVal, TIdx> A("examples/data/sparse_example.mtx", 1);
    DVector<TVal, TIdx> v{A.getCols(), 1.0};
    DVector<TVal, TIdx> u{A.getRows(), 1.0};

    GreedyVectorPartitioner<decltype(A), decltype(v)> pVecs(A, v, u);
    pVecs.partition();
    pVecs.localizeMatrix();

    u = A * v;
    JWLogVar(u);

    u = A * (v + v) + v;
    JWLogVar(u);

    // Dense matrices
    DMatrix<TVal, TIdx> D("examples/data/dense_example.mtx");
    DMatrix<TVal, TIdx> E("examples/data/dense_example.mtx");
    E.transpose();
    DMatrix<TVal, TIdx> F(D.getRows(), D.getRows());

    F = D * E;
    JWLogVar(D);
    JWLogVar(E);
    JWLogVar(F);
}
