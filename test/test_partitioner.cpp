#include <zee.h>
#include <cstdint.h>

using namespace Zee;

int main()
{
    // 1. Initialize a random n*m sparse matrix with set fillin
    int32_t n = 1000;
    int32_t m = 1000;
    double fill_in = 0.1;
    DSparseMatrix<double> A = rand(n, m, fill_in);

    // 2. Initialize a cyclic partitioner
    CyclicPartitioner<double> cycpart;

    // 3. Apply it to A
    cycpart.partition(A);

    // 4. Multiply A with some dense vector and store the result
    DDenseMatrix<double> v = rand(n, 1);
    DDenseMatrix<double> u;
    u = A * v;

    return 0;
}
