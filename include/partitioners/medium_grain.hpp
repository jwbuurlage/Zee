/* MG partitioner splits A into a 'column' and 'row' part, and then applies
 * a row partitioning. The identity blocks are implied, and used for computing
 * communication volume.
 *
 * There are four relevant matrices involved:
 *
 * A = A_r + A_c
 *
 * B = [ I   | A_r^T ]
 *     [-------------]
 *     [ A_c | I     ]
 *
 * We call B the "extended matrix".
 * This partitioner is a 'stateful iterative partitioner'.
 */

#pragma once

#include <zee.hpp>
using Zee::atomic_wrapper;

#include "kernighan_lin.hpp"

#include <random>

#include <vector>
using std::vector;

#include <memory>
using std::unique_ptr;
using std::make_unique;

#include <atomic>
using std::atomic;

template <class TMatrix = Zee::DSparseMatrix<double>>
class MGPartitioner : Zee::IterativePartitioner<TMatrix>
{
    public:
        using TIdx = typename TMatrix::index_type;
        using TVal = typename TMatrix::value_type;
        using TImage = typename TMatrix::image_type;
        using TTriplet = Zee::Triplet<TVal, TIdx>;

        MGPartitioner() {
            this->_procs = 2;
            this->_procs_in = 1;
        }

        virtual void initialize(TMatrix& A) override
        {
            this->split(A);
            _initialized = true;
        }

        bool locallyOptimal() const {
            return _locallyOptimal;
        }

        /** Split matrix A into A_r + A_c using score function */
        // maybe also construct B here
        void split(TMatrix& A)
        {
            // here we explicitely create a partitioning of A into
            // Ar and Ac
            //
            // We initialize a bitset for the nonzeros of A, and
            // let 0 be Ar and 1 be into Ac

            TIdx m = A.getRows();
            TIdx n = A.getCols();
            TIdx p = A.getProcs();

            _tripletInRow.resize(A.nonZeros());

            // we need to count the number of elements in each row
            // and column, use compute
            // NOTE: {row,col}_count_i = {row,col}_count[i % p][i / p];
            vector<vector<atomic_wrapper<TIdx>>> row_count(p,
                    vector<atomic_wrapper<TIdx>>(m / p + 1));
            vector<vector<atomic_wrapper<TIdx>>> col_count(p,
                    vector<atomic_wrapper<TIdx>>(n / p + 1));

            // FIXME: parallelize, using generalized compute
            TIdx s = 0;
            for (auto& pimg : A.getImages()) {
                for (auto key_count : pimg->getRowSet()) {
                    auto i = key_count.first;
                    auto num_r_i = key_count.second;
                    row_count[i % p][i / p].a += num_r_i;
                }

                for (auto key_count : pimg->getColSet()) {
                    auto j = key_count.first;
                    auto num_c_j = key_count.second;
                    col_count[j % p][j / p].a += num_c_j;
                }

                ++s;
            }

            // we now have distributed row/col counts.
            // time to split into bits
            TIdx cur = 0;
            for (auto& pimg : A.getImages()) {
                // for every nonzero in the image, we need to check if its
                // row count or col count is higher, and assign
                // the proper bool to the bitset
                for (auto& t : *pimg) {
                    auto i = t.row();
                    auto j = t.col();

                    auto s_r_i = row_count[i % p][i / p].a.load();
                    auto s_c_j = col_count[j % p][j / p].a.load();

                    // this is the score function
                    _tripletInRow[cur++] = (s_r_i < s_c_j);
                }

                ++s;
            }
        }

        void constructExtendedMatrix(TMatrix& A, TMatrix& B)
        {
            // ## Explicitely construct matrix B
            // we allow this to be distributed for generalizing to parallel
            // for now just put A^c and A^r on differnt procs (and diagonal
            // elements)
            // ALTERNATIVE: do column partitioning implicitely),
            // atm I think this is worse for a number of reasons:
            // - less flexibility (can easily switch (1D) partitioners for B
            //   if we explicitely construct it)
            // - cluttered code here
            std::vector<TTriplet> coefficients;

            auto& images = A.getMutableImages();

            std::set<TIdx> rowset;
            std::set<TIdx> colset;

            TIdx cur = 0;
            for (auto& pimg : images) {
                for (auto t : *pimg) {
                    bool tripletInAr = _tripletInRow[cur++];
                    auto bTriplet = TTriplet(
                                bRow(tripletInAr, t.row(), t.col(), A),
                                bCol(tripletInAr, t.row(), t.col(), A),
                                t.value());
                    rowset.insert(bTriplet.row());
                    colset.insert(bTriplet.col());
                    coefficients.push_back(bTriplet);
                }
            }

            for (TIdx i = 0; i < A.getRows() + A.getCols(); ++i) {
                if (colset.find(i) != colset.end() &&
                        rowset.find(i) != rowset.end()) {
                    coefficients.push_back(TTriplet(
                                i, i, (TVal)1));
                }
            }

            B.setFromTriplets(coefficients.begin(), coefficients.end());
        }

        inline TIdx bRow(bool tripletInAr,
                TIdx row,
                TIdx col,
                const TMatrix& A) const
        {
            if (tripletInAr)
                return col;
            else
                return A.getRows() + row;
        }

        inline TIdx bCol(bool tripletInAr,
                TIdx row,
                TIdx col,
                const TMatrix& A) const
        {
            if (tripletInAr)
                return A.getCols() + row;
            else
                return col;
        }

        /** The partitioning of the extended matrix B is used to repartition
         *  the matrix A */
        void inducePartitioning(TMatrix& A, TMatrix& B)
        {
            // we bipartition
            TIdx p = 2;

            vector<vector<int>> procForCol(p,
                    vector<int>(B.getRows() / p + 1, -1));
            TIdx s = 0;
            for (auto& subMatrix : B.getImages()) {
                for (auto& pColCount : subMatrix->getColSet()) {
                    auto col = pColCount.first;
                    procForCol[col % p][col / p] = s;
                }
                ++s;
            }

            vector<unique_ptr<TImage>> aNewImages;
            for (TIdx i = 0; i < p; ++i)
                aNewImages.push_back(std::make_unique<TImage>());

            auto& images = A.getMutableImages();
            TIdx cur = 0;
            // now we modify A
            for (auto& pimg : images) {
                for (auto triplet : *pimg) {
                    auto targetProc = 0;
                    if (_tripletInRow[cur++]) {
                        auto colProc = triplet.row() + A.getCols();
                        targetProc =
                            procForCol[colProc % p][colProc / p];
                    } else {
                        targetProc =
                            procForCol[triplet.col() % p][triplet.col() / p];
                    }
                    aNewImages[targetProc]->pushTriplet(triplet);
                }
            }

            A.resetImages(aNewImages);
        }

        /** For MG the initial partition relies on decomposition into A^c + A^r
         * using information about the distribution of A */
        TMatrix& partition(TMatrix& A) override
        {
            // FIXME: many of these operations can be done in place

            if (_initialized) {
                ZeeLogWarning << "Already applied an initial partitioning"
                    " instead refining the partitioning on A" << endLog;
                this->refine(A);
                return A;
            }

            if (!A.isInitialized()) {
                ZeeLogError << "MG: Trying to partition uninitialized matrix." << endLog;
                return A;
            }

            // PHASE 1: split A in two, and construct extended matrix B
            this->initialize(A);

            auto B = TMatrix(2 * A.getRows(), 2 * A.getCols());
            constructExtendedMatrix(A, B);

            // PHASE 2: call column partitioner
            // (we now partition B)
            // here we use Kernighan-Lin
            KernighanLin<decltype(B)> kerLin(B);
            kerLin.run();

            // PHASE 3: convert back partitioning to A
            // note that each image in B has a map proc -> col(s)
            // (distributed by proc) we want to construct from this a
            // (distributed) map col -> proc of size O(2n / p) per proc
            // lets just do this explicitely here
            inducePartitioning(A, B);

            return A;
        }


        /** IR the (bi-)partitioning on A */
        TMatrix& refine(TMatrix& A) override
        {
            // We only support IR on bi-partitionings
            if (A.getProcs() != 2) {
                ZeeLogError << "For now MG-IR only supports bipartitionings, "
                    "the matrix A has a " << A.getProcs() << "-way partitioning"
                    << endLog;
                return A;
            }

            if (!A.isInitialized()) {
                ZeeLogError << "MG: Trying to refine uninitialized matrix." << endLog;
                return A;
            }

            // We store the current communication volume such that we can
            // see if it gets improved
            auto priorVolume = A.communicationVolume();
            auto& aImages = A.getImages();

            auto B = TMatrix(2 * A.getRows(), 2 * A.getCols());

            std::vector<TTriplet> coefficients;
            coefficients.reserve((int)(A.nonZeros() + A.getRows() + A.getCols()));

            std::set<TIdx> rowset;
            std::set<TIdx> colset;

            std::fill(_tripletInRow.begin(), _tripletInRow.end(), false);

            // PHASE 1: We construct the matrix B according to the partitioning of A
            TIdx cur = 0;
            for (auto& t : *aImages[0]) {
                auto bTriplet = TTriplet(
                            bRow(_phaseRow, t.row(), t.col(), A),
                            bCol(_phaseRow, t.row(), t.col(), A),
                            t.value());
                if (_phaseRow) {
                    _tripletInRow[cur] = true;
                }
                rowset.insert(bTriplet.row());
                colset.insert(bTriplet.col());
                coefficients.push_back(bTriplet);
                cur++;
            }

            for (auto& t : *aImages[1]) {
                auto bTriplet = TTriplet(
                            bRow(!_phaseRow, t.row(), t.col(), A),
                            bCol(!_phaseRow, t.row(), t.col(), A),
                            t.value());
                rowset.insert(bTriplet.row());
                colset.insert(bTriplet.col());
                if (!_phaseRow) {
                    _tripletInRow[cur] = true;
                }
                coefficients.push_back(bTriplet);
                cur++;
            }

            for (TIdx i = 0; i < A.getCols() + A.getRows(); ++i) {
                if (colset.find(i) != colset.end() &&
                        rowset.find(i) != rowset.end()) {
                    coefficients.push_back(TTriplet(i, i, (TVal)1));
                }
            }

            // We want B to be bipartitioned.
            B.setDistributionScheme(Zee::partitioning_scheme::custom, 2);
            B.setDistributionFunction([&A] (TIdx row, TIdx col) {
                        if (col < A.getRows()) {
                            return 0;
                        }
                        return 1;
                    });
            B.setFromTriplets(coefficients.begin(), coefficients.end());

            if (B.isInitialized()) {
                // PHASE 2: Apply KL
                // FIXME: naming is off here, why B member? only used in initialize
                // FIXME: want to keep B inbetween runs, since it does things?
                KernighanLin<decltype(B)> kerLin(B);
                kerLin.run();

                // PHASE 3: retrieve partitioning for A
                inducePartitioning(A, B);
            }
            else {
                ZeeLogError << "Could not construct the extended matrix B."
                    << endLog;
            }

            if (A.communicationVolume() == priorVolume) {
                if (!_phaseRow && _rowOptimal) {
                    _locallyOptimal = true;
                } else if (_phaseRow) {
                    _rowOptimal = true;
                }
                _phaseRow = !_phaseRow;
            }
            else {
                if (!_phaseRow) {
                    _rowOptimal = false;
                }
            }

            return A;
        }

    private:
        vector<bool> _tripletInRow;
        bool _initialized = false;
        bool _locallyOptimal = false;
        bool _rowOptimal = false;
        bool _phaseRow = true;
};
