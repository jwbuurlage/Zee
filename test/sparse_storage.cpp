#include "catch.hpp"

#include "zee.hpp"

using TVal = Zee::default_scalar_type;
using TIdx = Zee::default_index_type;

TEST_CASE("moving non-zeros in storage", "[sparse storage]") {
    TIdx procs = 4;
    Zee::DSparseMatrix<> matrix(10, 11, 4);

    std::vector<Zee::Triplet<>> triplets = {{1, 1, 1}, {2, 2, 2}, {3, 3, 3},
                                            {4, 4, 4}, {5, 5, 5}, {6, 6, 6},
                                            {7, 7, 7}, {8, 8, 8}, {9, 9, 9}};

    TIdx k = 0;
    for (auto triplet : triplets) {
        matrix.pushTriplet((k++ % procs), triplet);
    }

    matrix.moveNonZero(0, 0, 1);
    REQUIRE(matrix[0].nonZeros() == 2);
    REQUIRE(matrix[1].nonZeros() == 3);
}
