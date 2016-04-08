#include "catch.hpp"

#include "zee.hpp"

TEST_CASE("explicitely set elements", "[loading matrices]") {
    Zee::DSparseMatrix<> matrix(10, 11, 1);
    matrix.pushTriplet(0, {1, 1, 1});
    matrix.pushTriplet(0, {2, 2, 2});
    REQUIRE(matrix.getCols() == 11);
    REQUIRE(matrix.getRows() == 10);
    REQUIRE(matrix.nonZeros() == 2);
}

TEST_CASE("special matrices", "[loading matrices]") {
    SECTION("load identity matrix") {
        unsigned int n = 10;
        unsigned int procs = 1;
        auto id = Zee::eye(n, procs);

        using TVal = decltype(id)::value_type;
        using TIdx = decltype(id)::index_type;

        Zee::DVector<TVal, TIdx> x(n, 0);
        Zee::DVector<TVal, TIdx> y(n, 0);
        Zee::GreedyVectorPartitioner<decltype(id), decltype(x)> part_vector(
            id, x, y);
        part_vector.partition();
        part_vector.localizeMatrix();

        y = id * x;
        REQUIRE(x == y);
    }

    SECTION("load random matrix") {
        unsigned int size = 100;
        unsigned int procs = 1;
        double density = 0.5;

        auto rand_matrix = Zee::rand(size, size, procs, density);
        REQUIRE(rand_matrix.nonZeros() > 0);
    }

}

TEST_CASE("dense matrices", "[loading matrices]") {
    REQUIRE(0 == 0);
}

TEST_CASE("matrix market", "[loading matrices]") {
    REQUIRE(0 == 0);
}

TEST_CASE("loading matrices for p > 1", "[loading matrices]") {
    REQUIRE(0 == 0);
}
