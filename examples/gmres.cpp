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
        TIdx maxIterations,
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

    DVector<TVal, TIdx> r = b;
    TVector e{center, A.getRows()};

    // We store V as a (centralized) pseudo matrix
    // I would think this is fine as long as m is small
    // This contains the orthogonal basis of our Krylov subspace
    std::vector<TVector> V(m, TVector(center, r.size()));

    // We store the upper Hessenberg H matrix containing increasingly
    // large vectors
    std::vector<TVector> H;
    H.reserve(m);
    for (auto i = 0; i < m; ++i) {
        H.push_back(TVector(center, i + 2));
    }

    // Similarly we store the matrix R
    std::vector<TVector> R(m, TVector(center, r.size()));

    // Store \hat{b}
    TVector bHat(center, A.getRows());

    // Additional variables used for the algorithm
    std::vector<TVal> c(m);
    std::vector<TVal> s(m);
    std::vector<TVal> y(m);

    std::vector<TVal> rhos;

    auto finished = false;
    auto nr = -1;
    for (auto run = 0; run < maxIterations; ++run) {
        // We construct the initial basis vector from the residual
        auto beta = r.norm();
        V[0] = r / beta;
        bHat.reset();
        bHat[0] = beta;

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
            for (int k = 0; k <= i; ++k) {
                H[i][k] = V[k].dot(w);
                w -= V[k] * H[i][k];
            }

            H[i][i + 1] = w.norm();

            if (i != m - 1)
                V[i + 1] = w / H[i][i + 1];

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
            c[i] = R[i][i] / delta;
            s[i] = H[i][i + 1] / delta;

            // update r_i and bHat
            R[i][i] = c[i] * R[i][i] + s[i] * H[i][i + 1];
            bHat[i + 1] = - s[i] * bHat[i];
            bHat[i] = c[i] * bHat[i];

            // update rho
            auto rho = std::abs(bHat[i + 1]);
            rhos.push_back(rho);

            bench.silence();

            // check if we are within tolerance level
            if (rho < tol) {
                nr = i;
                finished = true;
                break;
            }
        }

        // reconstruct x
        if (!finished) {
            nr = m - 1;
            y[nr] = bHat[nr] / R[nr][nr];
        }

        // subtract sum ( i dont know what i did )
        for (TIdx k = nr - 1; k >= 0; --k) {
            TVal sum = 0;
            for (TIdx i = k + 1; i <= nr; ++i) {
                sum += R[i][k] * y[i];
            }
            y[k] = (bHat[k] - sum) / R[k][k];
        }

        for (TIdx i = 0; i < nr; ++i) {
            x = x + y[i] * V[i];
        }

        r = b - A * x;

        if (finished) {
            break;
        }
    }

    // We plot the residuals
    auto p = Plotter<>();
    p["xlabel"] = "iterations";
    p["ylabel"] = "$\\rho$";
    p["yscale"] = "log";
    p["title"] = "GMRES: residual norm";
    p.addLine(rhos, "rhos");
    p.plot("residual_test", true);
}

} // namespace GMRES

int main()
{
    ZeeLogInfo << "-- Starting GMRES example" << endLog;

    auto bench = Benchmark("GMRES");
    bench.phase("initialize matrix");

    auto center = std::make_shared<Unpain::Center<int>>(4);

    // We initialize the matrix with cyclic distribution
    std::string matrix = "steam3";
    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx", 4);

    auto b = DVector<double, int>{center, A.getRows(), 1.0};
    b = A * b;

    // initial x is the zero vector
    auto x = DVector<double, int>{center, A.getCols()};

    // Start GMRES
    bench.phase("GMRES");
    GMRES::solve<double, int>(A, // Matrix
            b,                   // RHS vector
            x,                   // resulting guess for x
            100,                  // outer iterations
            100,                  // inner iterations
            1e-6);               // tolerance level

    DVector<double, int> c{center, A.getRows()};
    DVector<double, int> r{center, A.getRows()};
    c = A * x;
    r = b - c;
    ZeeLogVar(r.norm() / b.norm());

    return 0;
}
