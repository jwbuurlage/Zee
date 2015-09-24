#include <zee.hpp>

using namespace Zee;

int main()
{
    int N = 4;
    int M = 16;
    int core_block_size = 8;
    int matrix_size = core_block_size * N * M;

    // Dense matrices
    DMatrix<> A(matrix_size, matrix_size);
    DMatrix<> B(matrix_size, matrix_size);

    for (int i = 0; i < matrix_size; ++i)
        for (int j = 0; j < matrix_size; ++j) {
            A.at(i, j) = (float)i / 10.0f;
            B.at(i, j) = (float)j / 10.0f;
        }
    A.spy("A");
    B.spy("B");

    DMatrix<> C(A.getRows(), A.getRows());

    C = A * B;
    ZeeLogVar(matrix_size);
    C.spy("C");

    return 0;
}
