#include <sstream>

#include "catch.hpp"

#include "zee.hpp"

using TVal = Zee::default_scalar_type;
using TIdx = Zee::default_index_type;

TEST_CASE("explicitely set elements", "[loading matrices]") {
    Zee::DSparseMatrix<> matrix(10, 11, 1);

    std::vector<Zee::Triplet<>> triplets = {
        {1, 1, 1}, {2, 2, 2},
    };

    for (auto triplet : triplets) {
        matrix.pushTriplet(0, triplet);
    }

    REQUIRE(matrix.getCols() == 11);
    REQUIRE(matrix.getRows() == 10);
    REQUIRE(matrix.nonZeros() == 2);
}

TEST_CASE("special matrices", "[loading matrices]") {
    SECTION("load identity matrix") {
        unsigned int n = 10;
        unsigned int procs = 1;
        auto id = Zee::eye(n, procs);

        using Val = decltype(id)::value_type;
        using Idx = decltype(id)::index_type;

        Zee::DVector<Val, Idx> x(n, 0);
        Zee::DVector<Val, Idx> y(n, 0);
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
    SECTION("set directly") {
        Zee::DMatrix<> matrix(10, 11);
        matrix.at(1, 1) = 1;
        matrix.at(6, 5) = 6;
        REQUIRE(matrix.at(1, 1) == 1);
        REQUIRE(matrix.at(6, 5) == 6);
        REQUIRE(matrix.at(0, 0) == 0);
    }

    SECTION("set from values") {
        Zee::DMatrix<> matrix(4, 3);
        std::vector<TVal> values = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
        matrix.setFromValues(values.begin(), values.end());
        REQUIRE(matrix.at(2, 1) == 8);
    }
}

TEST_CASE("matrix market", "[loading matrices]") {
    SECTION("sparse") {
        SECTION("basic loading") {
            Zee::DSparseMatrix<> matrix{"test/mtx/sparse_example.mtx", 1};
            REQUIRE(matrix.nonZeros() == 3);

            Zee::DSparseMatrix<> square_matrix{
                "test/mtx/square_sparse_example.mtx", 1};
            REQUIRE(square_matrix.nonZeros() == 4);
        }

        SECTION("symmetric") {
            Zee::DSparseMatrix<> symmetric_matrix{
                "test/mtx/symmetric_sparse_example.mtx", 1};
            REQUIRE(symmetric_matrix.nonZeros() == 5);
        }

        SECTION("skew-symmetric") {
            Zee::DSparseMatrix<> skew_symmetric_matrix{
                "test/mtx/skew_symmetric_sparse_example.mtx", 1};
            REQUIRE(skew_symmetric_matrix.nonZeros() == 5);
        }

        SECTION("save-and-load") {
            // FIXME: require that saving and loading matrix market gives the
            // same matrix (up to sparsity)

            // FIXME dont want this to depend on file system, add spy-as-string
            Zee::DSparseMatrix<> matrix{"test/mtx/sparse_example.mtx"};
            std::stringstream ss;
            matrix.spyToStream(ss);

            Zee::DSparseMatrix<> matrixSaveLoad;
            Zee::matrix_market::loadFromStream(ss, matrixSaveLoad);

            REQUIRE(matrix == matrixSaveLoad);
        }
    }
    SECTION("dense") {
        Zee::DMatrix<TVal, TIdx> matrix{"test/mtx/dense_example.mtx"};
        REQUIRE(matrix.at(2, 0) == 3);
    }
}

TEST_CASE("loading matrices for p > 1", "[loading matrices]") {
    TIdx procs = 4;
    Zee::DSparseMatrix<> matrix(10, 11, 4);

    std::vector<Zee::Triplet<>> triplets = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3},
                                            {4, 4, 4}, {5, 5, 5}, {6, 6, 6},
                                            {7, 7, 7}, {8, 8, 8}, {9, 9, 9}};

    TIdx k = 0;
    for (auto triplet : triplets) {
        matrix.pushTriplet((k++ % procs), triplet);
    }

    REQUIRE(matrix.getCols() == 11);
    REQUIRE(matrix.getRows() == 10);
    REQUIRE(matrix[0].nonZeros() == 3);
    REQUIRE(matrix[1].nonZeros() == 2);
    REQUIRE(matrix[2].nonZeros() == 2);
    REQUIRE(matrix[3].nonZeros() == 2);
}
