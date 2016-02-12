#include <string>
#include <cmath>

#include <zee.hpp>

using namespace Zee;

int main(int argc, char* argv[]) {
    using TVal = double;
    using TIdx = uint32_t;

    auto args = ArgParse();
    args.addOption("--matrices", "list of matrices to partition", true);
    if (!args.parse(argc, argv))
        return -1;

    auto matrices = args.asList("--matrices");

    // We initialize the matrix with cyclic distribution
    for (auto matrix : matrices) {
        auto A = DSparseMatrix<TVal, TIdx>("data/matrices/" + matrix + ".mtx", 1);

        auto b = DVector<TVal, TIdx>{A.getRows(), 1.0};
        GreedyVectorPartitioner<decltype(A), decltype(b)> pVecs(A, b, b);
        pVecs.partition();
        pVecs.localizeMatrix();

        b = A * b;

        // initial x is the zero vector
        auto x = DVector<TVal, TIdx>{A.getCols()};

        // Start GMRES
        GMRES::solve<TVal, TIdx>(A, // Matrix
                b,                  // RHS vector
                x,                  // resulting guess for x
                1,                // outer iterations
                50,                // inner iterations
                1e-6,               // tolerance level
                true);             // plot residuals

        DVector<TVal, TIdx> c{A.getRows()};
        DVector<TVal, TIdx> r{A.getRows()};
        c = A * x;
        r = b - c;
        ZeeLogVar(r.norm() / b.norm());
    }

    return 0;
}
