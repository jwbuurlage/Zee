#include <zee.hpp>

using namespace Zee;

int main()
{

    int matrix_size = 16;

    // Dense matrices
    DMatrix<> A(matrix_size, matrix_size);
    DMatrix<> B(matrix_size, matrix_size);

    for (int i = 0; i < matrix_size; ++i)
        for (int j = 0; j < matrix_size; ++j) {
            A.at(i, j) = (float)i;
            B.at(i, j) = (float)j;
        }
    A.spy("A");
    B.spy("B");

    DMatrix<> C(A.getRows(), A.getRows());

    C = A * B;
    ZeeLogVar(matrix_size);
    ZeeLogVar(C);
    C.spy("C");

    return 0;
}
