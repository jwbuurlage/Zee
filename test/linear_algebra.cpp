#include "catch.hpp"

#include <zee.hpp>

using TIdx = uint32_t;
using TVal = float;

TIdx size = 4;

Zee::DVector<TVal, TIdx> x{size, 1.0};
Zee::DVector<TVal, TIdx> y{size, 1.0};
Zee::DVector<TVal, TIdx> z{size, 1.0};

TEST_CASE("vector operations", "[linear algebra]") {
    REQUIRE(x.size() == size);

    SECTION("we can add vectors") {
        z = x + y;
        REQUIRE(z[0] == 2.0f);
    }

    SECTION("we can subtract vectors") {
        z = x - y;
        REQUIRE(z[0] == 0.0f);
    }

    SECTION("we can add multiple vectors in various ways") {
        z = x + x + x;
        REQUIRE(z[0] == 3.0f);

        z = (x + x) + x;
        REQUIRE(z[0] == 3.0f);

        z = x + (x + x);
        REQUIRE(z[0] == 3.0f);

        z = (x + x) + (x + x);
        REQUIRE(z[0] == 4.0f);
    }

    SECTION("we can perform scalar multiplication on vectors") {
        z = 2.5f * x;
        REQUIRE(z[0] == 2.5f);

        z = x * 2.5f;
        REQUIRE(z[0] == 2.5f);

        z = x / 2.0f;
        REQUIRE(z[0] == 0.5f);
    }

    SECTION("we can mix vector operations") {
        z = 2.5f * (x + y) / 0.5 + x + x + (x - y);
        REQUIRE(z[0] == 12.0f);
    }
}

Zee::DSparseMatrix<TVal, TIdx> A{"mtx/sparse_example.mtx", 1};
Zee::DMatrix<TVal, TIdx> B{"mtx/dense_example.mtx"};
Zee::DMatrix<TVal, TIdx> C = B;
Zee::DMatrix<TVal, TIdx> D{4, 4};

TEST_CASE("matrix operations", "[linear algebra]") {
    SECTION("we can load matrix market format") {
        REQUIRE(A.getRows() == 3);
        REQUIRE(A.getCols() == 4);

        REQUIRE(B.getRows() == 3);
        REQUIRE(B.getCols() == 4);

        REQUIRE(B.at(2, 3) == 12.0f);
    }

    SECTION("we can perform a spmv") {
        z = A * x;
        REQUIRE(z[1] == 2.0f);
        REQUIRE(z.size() == 3);
    }

    SECTION("we can mix spmv and vectors") {
        z = A * (x + x) * 2.0f;
        REQUIRE(z[1] == 8.0f);
        REQUIRE(z.size() == 3);
    }

    SECTION("dense operations") {
        SECTION("we can transpose (dense) matrices") {
            C.transpose();
            REQUIRE(C.getRows() == B.getCols());
            REQUIRE(C.getCols() == B.getRows());
            REQUIRE(C.at(3, 2) == 12.0f);
        }

        SECTION("we can multiply dense matrices") {
            D = B * C;
            REQUIRE(D.at(2, 2) == 270.0f);
        }
    }
}
