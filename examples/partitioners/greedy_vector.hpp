#include <zee.hpp>

template <class TMatrix = DSparseMatrix<>, class TVector = DVector<>>
class GreedyVectorPartitioner : Zee::VectorPartitioner<TMatrix, TVector>
{
    using TIdx = TMatrix::index_type;

    public:
        GreedyVectorPartitioner(TMatrix& A, TVector& v, TVector& u)
            : Zee::VectorPartitioner(A, v, u)
        {} 

        void partition() override {
            // We assign P_v(k) = P_u(k) according to the following scheme:
            // 1. Assign to P_A(k, k)
            // 2. If empty, assign to target to some proc in intersection of
            //    P(k, *) and P(*, k), greedily balancing the assignment load
            // 3. If empty, assign to target to some proc in union of P(k, *)
            //    and P(*, k), greedily balancing the assignment load

            // We assume the matrix A_ is partitioned
            ZeeAssert(A.isInitialized());

            // need to assert matrix A_ is square for this particular solver
            // TODO figure out how does CGLS work, and how to generalize
            // dist(u) = dist(v) for those solvers
            ZeeAssert(A.getRows() == A.getCols());
            ZeeAssert(v.size() == A.getRows());
            ZeeAssert(u.size() == v.size());

            int n = v.size();
            int p = A.getProcs();

            // We need to know if diagonal is nonzero
            // request diagonals from images
            // can preallocate (distributed) array of size n
            // let them write into that.. or use message queue system
            vector<int> diagonalTargets{-1, n};
            // assume diagonal locations are known
            // (in diagonal targets)

            for (TIdx i = 0; i < v.size(); ++i) {
                // FIXME implement vector.reassign
                if (diagonalTargets[i] != -1) {
                    v.reassign(i, diagonalTargets[i]);
                    u.reassign(i, diagonalTargets[i]);
                }
                else if (/* intersection != empty */) {
                    // compute intersection
                    // ask processors if they have non-empty row/column k
                    // excluding diagonal
                    // greedily balance
                    // hard
                }
                else {
                    // compute union
                    // greedily balance
                    // hard
                }
            }
        };
};
