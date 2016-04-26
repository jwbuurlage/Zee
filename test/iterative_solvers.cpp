#include "catch.hpp"

#include "zee.hpp"

TEST_CASE("we can solve a non-symmetric system", "[solvers]") {
    Zee::DSparseMatrix<> matrix{"test/mtx/ex24.mtx", 1};

    auto n = matrix.getCols();

    auto ones = Zee::DVector<>{n, 1.0};
    auto b = Zee::DVector<>{n, 0.0};

    Zee::GreedyVectorPartitioner<decltype(matrix), decltype(b)>
        vector_partitioner(matrix, ones, b);
    vector_partitioner.partition();
    vector_partitioner.localizeMatrix();

    b = matrix * ones;

    auto x = Zee::DVector<>{n, 0.0};

    auto tol = 1e-2;

    Zee::GMRES::solve<Zee::default_scalar_type, Zee::default_index_type>(
        matrix, b, x, 1, n, tol);

    auto r = Zee::DVector<>{n, 0.0};
    r = x - ones;

    REQUIRE(r.norm() / ones.norm() < 100 * tol);
}
