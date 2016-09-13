#include <algorithm>
#include <vector>

#include "vector_partitioner.hpp"

namespace Zee {

template <class TMatrix = DSparseMatrix<>, class TVector = DVector<>>
class GreedyVectorPartitioner : public VectorPartitioner<TMatrix, TVector> {
    using TIdx = typename TMatrix::index_type;

   public:
    GreedyVectorPartitioner(TMatrix& A, TVector& v, TVector& u)
        : VectorPartitioner<TMatrix, TVector>(A, v, u) {}

    void partition() override {
        // Aliases
        auto& A = this->A_;
        auto& u = this->u_;
        auto& v = this->v_;

        using TIdx = typename TMatrix::index_type;

        // We assign P_v(k) = P_u(k) according to the following scheme:
        // 1. Assign to P_A(k, k)
        // 2. If empty, assign to target to some proc in intersection of
        //    P(k, *) and P(*, k), greedily balancing the assignment load
        // 3. If empty, assign to target to some proc in union of P(k, *)
        //    and P(*, k), greedily balancing the assignment load

        // We assume the matrix A_ is partitioned
        JWAssert(A.isInitialized());

        // need to assert matrix A_ is square for this particular solver
        // TODO figure out how does CGLS work, and how to generalize
        // dist(u) = dist(v) for those solvers
        JWAssert(A.getRows() == A.getCols());
        JWAssert(v.size() == A.getRows());
        JWAssert(u.size() == v.size());

        TIdx n = v.size();
        TIdx p = A.getProcs();

        std::vector<std::set<TIdx>> processorsInRow(n);
        std::vector<std::set<TIdx>> processorsInCol(n);

        // We need to know if diagonal is nonzero
        // FIXME: To parallelize: request diagonals from images
        // can preallocate (distributed) array of size n
        // let them write into that.. or use message queue system
        std::vector<TIdx> diagonalTargets(n, p);

        TIdx s = 0;
        for (auto& img : A.getImages()) {
            for (auto& triplet : *img) {
                if (triplet.col() == triplet.row()) {
                    diagonalTargets[triplet.col()] = s;
                }
                processorsInRow[triplet.row()].insert(s);
                processorsInCol[triplet.col()].insert(s);
            }
            ++s;
        }

        std::vector<TIdx> elementCount(p, 0);

        for (TIdx i = 0; i < v.size(); ++i) {
            if (diagonalTargets[i] != p) {
                v.reassign(i, diagonalTargets[i]);
                u.reassign(i, diagonalTargets[i]);
                elementCount[diagonalTargets[i]]++;
            } else {
                std::set<TIdx> lookUpSet;  // intersection or union
                std::set_intersection(
                    processorsInRow[i].begin(), processorsInRow[i].end(),
                    processorsInCol[i].begin(), processorsInCol[i].end(),
                    std::inserter(lookUpSet, lookUpSet.begin()));
                if (lookUpSet.empty()) {
                    std::set_union(
                        processorsInRow[i].begin(), processorsInRow[i].end(),
                        processorsInCol[i].begin(), processorsInCol[i].end(),
                        std::inserter(lookUpSet, lookUpSet.begin()));
                }
                auto minElement = std::min_element(
                    lookUpSet.begin(), lookUpSet.end(),
                    [&, elementCount](const TIdx& a, const TIdx& b) {
                        return elementCount[a] < elementCount[b];
                    });

                TIdx lightestProc = *minElement;
                v.reassign(i, lightestProc);
                u.reassign(i, lightestProc);
                elementCount[lightestProc]++;
            }
        }
    };
};

}  // namespace Zee
