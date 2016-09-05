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

TEST_CASE("hyper-pulp gives a valid partitioning", "[partitioning]") {
    SECTION("matrix sparsity does not change") {
        Zee::DSparseMatrix<TVal, TIdx> original{"test/mtx/tomo_matrix.mtx", 1};
        Zee::DSparseMatrix<TVal, TIdx> matrix{"test/mtx/tomo_matrix.mtx", 1};

        // partition matrix, allow for initial imbalance
        Zee::PulpPartitioner<decltype(matrix)> pulp(matrix, 2, 2.0);
        pulp.initialize(matrix);
        pulp.initialPartitioning(matrix.getCols());
        auto com_vol = matrix.communicationVolume();
        while(true) {
            pulp.refineWithIterations(matrix.getCols());
            if (com_vol == matrix.communicationVolume())
                break;
            com_vol = matrix.communicationVolume();
        }

        CHECK(original == matrix);
    }
}
