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
        DVector<TVal, TIdx>& x)
{
    using TVector = DVector<TVal, TIdx>;

    const auto& center = A.getCenter();

    auto r = b;
    TVector e{center, A.getRows()};

    // FIXME for restarts.. not applicable yet
    //auto maxIterations = 10;

    // For convenience we call m the number of inner iterations,
    // which is also the maximum size of for example the H matrix
    auto m = 100;

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
    ZeeLogVar(beta);

    V.push_back(TVector(center, r.size()));
    V[0] = r / beta;

    // We construct \hat{b}
    TVector bHat(center, A.getRows());
    bHat[0] = beta;

    // Additional variables used for the algorithm
    std::vector<TVal> c;
    std::vector<TVal> s;
    std::vector<TVal> y;

    // We run for i [0, m)
    for (auto i = 0; i < m; ++i)
    {
        // We introduce a new basis vector which we will orthogonalize
        // using modified Gramm-Schmidt
        TVector w(center, A.getRows());
        w = A * V[i];

        // H[i]
        H.push_back(TVector(center, i + 2));
        for (int k = 0; k <= i; ++k) {
            // update h_ki
            H[i][k] = V[k].dot(w);
            w -= V[k] * H[i][k];
        }

        H[i][i + 1] = w.norm();

        // V_(i + 1)
        V.push_back(TVector(center, r.size()));
        V[i + 1] = w / H[i][i + 1];

        // R_i
        R.push_back(TVector(center, r.size()));
        R[i][0] = H[i][0];

        // Givens rotation
        for (int k = 1; k <= i; ++k) {
            auto gamma = c[k - 1] * R[i][k - 1] + s[k - 1] * H[i][k];
            R[i][k] = c[k - 1] * H[i][k] - s[k - 1] * R[i][k - 1];
            R[i][k - 1] = gamma;
        }

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
        ZeeLogInfo << "||b - A x" << i << "|| = " << rho << endLog;
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
    auto center = std::make_shared<Unpain::Center<int>>(2);
    std::string matrix = "karate";

    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx");
    auto b = DVector<double, int>{center, A.getCols(), 1.0};

    bench.phase("partition");

    // We partition the matrix with medium grain
    MGPartitioner<decltype(A)> partitioner;
    partitioner.partition(A);
    for (auto iter = 0; iter < 100; ++iter) {
        partitioner.refine(A);
    }
    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());

    // initial x is the zero vector
    auto x = DVector<double, int>{center, A.getRows()};

    bench.phase("GMRES");
    // Start GMRES
    GMRES::solve<double, int>(A, b, x);

//    auto r = b - A * x;
//    ZeeLogInfo << "|| b - A x || = : " << r.norm() << endLog;

    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());

    return 0;
}
