#include "catch.hpp"

#include "zee.hpp"

using TVal = Zee::default_scalar_type;
using TIdx = Zee::default_index_type;

TEST_CASE("cyclic partitioning", "[partitioning]") {
    TIdx procs = 2;

    // initialize the matrix from file
    Zee::DSparseMatrix<TVal, TIdx> matrix{"test/mtx/sparse_example.mtx", 1};
    Zee::CyclicPartitioner<decltype(matrix)> cyclic_part(procs, Zee::CyclicType::row);
    cyclic_part.partition(matrix);

    REQUIRE(matrix.getProcs() == procs);
}
