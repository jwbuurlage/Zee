#include <zee.hpp>

#include <vector>

namespace GMRES
{

template<typename TVal, typename TIdx>
void solve(Zee::DSparseMatrix<TVal, TIdx>& A,
        const Zee::DVector<TVal, TIdx>& b,
        Zee::DVector<TVal, TIdx>& x,
        TIdx maxIterations,
        TIdx m,
        TVal tol,
        bool plotResiduals = false,
        bool benchmark = false)
{
    ZeeLogInfo << "Solving Ax = b for system of size " << A.getRows() <<
        " x " << A.getCols() << " with " << A.nonZeros() << " non-zeros" << endLog;

    // FIXME pointless to make bench if we dont benchmark
    auto bench = Zee::Benchmark("GMRES");
    if (!benchmark)
        bench.silence();

    ZeeAssert(A.getRows() == b.size());

    using TVector = Zee::DVector<TVal, TIdx>;

    // make sure m is not larger than RHS vector
    m = std::min(A.getRows(), m);

    TVector r = b;
    TVector e{A.getRows()};

    // We store V as a (centralized) pseudo matrix
    // I would think this is fine as long as m is small
    // This contains the orthogonal basis of our Krylov subspace
    std::vector<TVector> V(m, TVector(r.size()));

    // We store the upper Hessenberg H matrix containing increasingly
    // large vectors
    std::vector<TVector> H;
    H.reserve(m);
    for (TIdx i = 0; i < m; ++i) {
        H.push_back(TVector(i + 2));
    }

    // Similarly we store the matrix R
    std::vector<TVector> R(m, TVector(r.size()));

    // Store \hat{b}
    TVector bHat(A.getRows());

    // Additional variables used for the algorithm
    std::vector<TVal> c(m);
    std::vector<TVal> s(m);
    std::vector<TVal> y(m);

    std::vector<TVal> rhos;

    auto finished = false;
    TIdx nr = 0;
    for (TIdx run = 0; run < maxIterations; ++run) {
        // We construct the initial basis vector from the residual
        auto beta = r.norm();
        V[0] = r / beta;
        bHat.reset();
        bHat[0] = beta;

        // We run for i [0, m)
        for (TIdx i = 0; i < m; ++i)
        {

            // We introduce a new basis vector which we will orthogonalize
            // using modified Gramm-Schmidt
            TVector w(A.getRows());

            w = A * V[i];

            for (TIdx k = 0; k <= i; ++k) {
                H[i][k] = V[k].dot(w);
                w -= V[k] * H[i][k];
            }

            H[i][i + 1] = w.norm();

            if (i != m - 1)
                V[i + 1] = w / H[i][i + 1];

            R[i][0] = H[i][0];

            // Givens rotation
            for (TIdx k = 1; k <= i; ++k) {
                auto gamma = c[k - 1] * R[i][k - 1] + s[k - 1] * H[i][k];
                R[i][k] = c[k - 1] * H[i][k] - s[k - 1] * R[i][k - 1];
                R[i][k - 1] = gamma;
            }

            // update c and s
            auto delta = sqrt(R[i][i] * R[i][i] + H[i][i + 1] * H[i][i + 1]);

            if (delta < 1e-6) {
                nr = i;
                finished = true;
                break;
            }

            c[i] = R[i][i] / delta;
            s[i] = H[i][i + 1] / delta;

            // update r_i and bHat
            R[i][i] = c[i] * R[i][i] + s[i] * H[i][i + 1];
            bHat[i + 1] = - s[i] * bHat[i];
            bHat[i] = c[i] * bHat[i];

            // update rho
            auto rho = std::abs(bHat[i + 1]);
            rhos.push_back(rho);
            ZeeLogVar(rho);

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

        // subtract sum
        for (TIdx k = nr - 1;; --k) {
            TVal sum = 0;
            for (TIdx i = k + 1; i <= nr; ++i) {
                sum += R[i][k] * y[i];
            }

            y[k] = (bHat[k] - sum) / R[k][k];

            if (k == 0)
                break;
        }

        for (TIdx i = 0; i < nr; ++i) {
            x = x + y[i] * V[i];
        }

        r = b - A * x;

        if (finished) {
            break;
        }
    }

    // we finalize the benchmark
    if (benchmark)
        bench.finish();

    if (plotResiduals) {
        ZeeLogVar(rhos);
        // We plot the residuals
        auto p = Zee::Plotter<>();
        p["xlabel"] = "iterations";
        p["ylabel"] = "$\\rho$";
        p["yscale"] = "log";
        p["title"] = "GMRES: residual norm";
        p.addLine(rhos, "rhos");
        p.plot("residual_test", true);
    }
}

} // namespace GMRES
