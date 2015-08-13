#include <string>
#include <cmath>

#include <unpain_cpp.hpp>
#include <zee.hpp>

#include "medium_grain/medium_grain.hpp"

using namespace Zee;

namespace GMRES
{

template<typename TVal, typename TIdx>
void solve(const DSparseMatrix<TVal, TIdx>& A,
        const DVector<TVal, TIdx>& b,
        DVector<TVal, TIdx>& x,
        TIdx m,
        TVal tol)
{
    ZeeLogInfo << "Solving Ax = b for system of size " << A.getRows() <<
        " x " << A.getCols() << " with " << A.nonZeros() << " non-zeros" << endLog;

    ZeeAssert(A.getRows() == b.size());

    using TVector = DVector<TVal, TIdx>;

    // make sure m is not larger than RHS vector
    m = std::min(A.getRows(), m);

    const auto& center = A.getCenter();

    auto r = b;
    TVector e{center, A.getRows()};

    // FIXME support restarts..
    auto maxIterations = 1;

    // We store V as a (centralized) pseudo matrix
    // I would think this is fine as long as m is small
    // This contains the orthogonal basis of our Krylov subspace
    std::vector<TVector> V;

    // We store the upper Hessenberg H matrix containing increasingly
    // large vectors
    std::vector<TVector> H;

    // Similarly we store the matrix R
    std::vector<TVector> R;

    // We construct the initial basis vector from the residual
    auto beta = r.norm();

    V.push_back(TVector(center, r.size()));
    V[0] = r / beta;

    // We construct \hat{b}
    TVector bHat(center, A.getRows());
    bHat[0] = beta;

    // Additional variables used for the algorithm
    std::vector<TVal> c;
    std::vector<TVal> s;
    std::vector<TVal> y;

    for (auto run = 0; run < maxIterations; ++run) {
        // We run for i [0, m)
        for (auto i = 0; i < m; ++i)
        {
            auto bench = Benchmark("GMRES inner loop");

            bench.phase("SpMV");

            // We introduce a new basis vector which we will orthogonalize
            // using modified Gramm-Schmidt
            TVector w(center, A.getRows());
            w = A * V[i];

            bench.phase("Gram-Schmidt");

            H.push_back(TVector(center, i + 2));
            for (int k = 0; k <= i; ++k) {
                // update h_ki
                H[i][k] = V[k].dot(w);
                w -= V[k] * H[i][k];
            }

            H[i][i + 1] = w.norm();

            bench.phase("V, R");

            V.push_back(TVector(center, r.size()));
            V[i + 1] = w / H[i][i + 1];

            R.push_back(TVector(center, r.size()));
            R[i][0] = H[i][0];

            bench.phase("Givens rotation");
            // Givens rotation
            for (int k = 1; k <= i; ++k) {
                auto gamma = c[k - 1] * R[i][k - 1] + s[k - 1] * H[i][k];
                R[i][k] = c[k - 1] * H[i][k] - s[k - 1] * R[i][k - 1];
                R[i][k - 1] = gamma;
            }

            bench.phase("Update residual");
            // update c and s
            auto delta = sqrt(R[i][i] * R[i][i] + H[i][i + 1] * H[i][i + 1]);
            c.push_back(R[i][i] / delta);
            s.push_back(H[i][i + 1] / delta);

            // update r_i and bHat
            R[i][i] = c[i] * R[i][i] + s[i] * H[i][i + 1];
            bHat[i + 1] = - s[i] * bHat[i];
            bHat[i] = c[i] * bHat[i];

            auto rho = std::abs(bHat[i + 1]);

            // rho is equal to the current error
            //ZeeLogInfo << "||b - A x" << i << "|| = " << rho << endLog;

            bench.silence();

            if (rho < tol)
                break;
        }
    }

    // TODO reconstruct x
}

} // namespace GMRES

int main()
{
    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    auto bench = Benchmark("GMRES");

    bench.phase("initialize matrix");

    // We initialize the matrix and rhs vector
    auto center = std::make_shared<Unpain::Center<int>>(4);
    std::string matrix = "fpga_dcop_05";

    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx", 1);
    A.spy();

    auto b = DVector<double, int>{center, A.getRows(), 1.0};
    b = A * b;

    bench.phase("partition");

    // We partition the matrix with medium grain
    //MGPartitioner<decltype(A)> partitioner;
    //partitioner.partition(A);
    //while (!partitioner.locallyOptimal()) {
    //    partitioner.refine(A);
    //}
    //ZeeLogVar(A.communicationVolume());
    //ZeeLogVar(A.loadImbalance());

    // initial x is the zero vector
    auto x = DVector<double, int>{center, A.getCols()};

    // Start GMRES
    bench.phase("GMRES");
    GMRES::solve<double, int>(A,        // Matrix
            b,                          // RHS vector
            x,                          // resulting guess for x
            50,                          // inner iterations
            1e-6);                      // tolerance level

    return 0;
}
