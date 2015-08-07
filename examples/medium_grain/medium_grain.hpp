/* MG partitioner splits A into a 'column' and 'row' part, and then applies
 * a row partitioning. The identity blocks are implied, and used for computing
 * communication volume.
 *
 * A = A_r + A_c
 * B = [ I & A_c^T \\ A_r & I ]
 *
 * We call B the "extended matrix".
 * This partitioner is a stateful iterative partitioner.
 */

#include <zee.hpp>
using Zee::atomic_wrapper;

#include "../kernighan_lin/kernighan_lin.hpp"

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

        MGPartitioner() {
            this->_procs = 2;
            this->_procs_in = 1;
        }

        virtual void initialize(TMatrix& A) override
        {
            this->split(A);
            initialized = true;
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

            auto m = A.rows();
            auto n = A.cols();
            auto p = A.procs();

            _bitInRow = make_unique<vector<vector<atomic_wrapper<bool>>>>(p);
            for (TIdx s = 0; s < p; ++s)
                (*_bitInRow)[s].resize((*A.getImages()[s]).nonZeros());

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
                    row_count[i % p][i / p]._a += num_r_i;
                }

                for (auto key_count : pimg->getColSet()) {
                    auto j = key_count.first;
                    auto num_c_j = key_count.second;
                    col_count[j % p][j / p]._a += num_c_j;
                }

                ++s;
            }

            // we now have distributed row/col counts.
            // time to split into bits
            s = 0;
            for (auto& pimg : A.getImages()) {
                // for every nonzero in the image, we need to check if its
                // row count or col count is higher, and assign
                // the proper bool to the bitset
                TIdx cur = 0;
                for (auto& t : *pimg) {
                    auto i = t.row();
                    auto j = t.col();

                    auto s_r_i = row_count[i % p][i / p]._a.load();
                    auto s_c_j = col_count[j % p][j / p]._a.load();

                    // this is the score function
                    (*_bitInRow)[s][cur++]._a = (s_r_i < s_c_j);
                }

                ++s;
            }
        }

        void constructExtendedMatrix(TMatrix& A, TMatrix& B)
        {
            // we bipartition
            auto p = 2;

            // ## Explicitely construct matrix B
            // we allow this to be distributed for generalizing to parallel
            // for now just put A^c and A^r on differnt procs (and diagonal
            // elements)
            // ALTERNATIVE: do column partitioning implicitely),
            // atm I think this is worse for a number of reasons:
            // - less flexibility (can easily switch (1D) partitioners for B
            //   if we explicitely construct it)
            // - cluttered code here
            vector<unique_ptr<TImage>> new_images;
            for (TIdx i = 0; i < p; ++i)
                new_images.push_back(std::make_unique<TImage>());

            auto& images = A.getMutableImages();

            TIdx s = 0;
            for (auto& pimg : images) {
                TIdx cur = 0;
                for (auto t : *pimg) {
                    if ((*_bitInRow)[s][cur++]._a)
                        new_images[0]->pushTriplet(
                                Zee::Triplet<TVal, TIdx>(
                                    t.col(), A.cols() + t.row(), t.value()));
                    else
                        new_images[1]->pushTriplet(
                                Zee::Triplet<TVal, TIdx>(
                                    A.rows() + t.row(), t.col(), t.value()));
                }
                ++s;
            }

            for (TIdx i = 0; i < A.rows() + A.cols(); ++i) {
                if (i < A.cols()) {
                    new_images[0]->pushTriplet(
                            Zee::Triplet<TVal, TIdx>(
                                i, i, (TVal)1));
                } else {
                    new_images[1]->pushTriplet(
                            Zee::Triplet<TVal, TIdx>(
                                i, i, (TVal)1));
                }
            }

            B.resetImages(new_images);
        }

        inline TIdx bRow(bool tripletInAr,
                TIdx row,
                TIdx col,
                const TMatrix& A) const
        {
            if (tripletInAr)
                return col;
            else
                return A.rows() + row;
        }

        inline TIdx bCol(bool tripletInAr,
                TIdx row,
                TIdx col,
                const TMatrix& A) const
        {
            if (tripletInAr)
                return A.cols() + row;
            else
                return col;
        }

        /** The partitioning of the extended matrix B is used to repartition
         *  the matrix A */
        void inducePartitioning(TMatrix& A, TMatrix& B)
        {
            // we bipartition
            auto p = 2;

            vector<vector<int>> procForCol(p,
                    vector<int>(B.rows() / p + 1, -1));
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
            // now we modify A
            s = 0;
            for (auto& pimg : images) {
                TIdx cur = 0;
                for (auto triplet : *pimg) {
                    auto targetProc = 0;
                    if ((*_bitInRow)[s][cur++]._a) {
                        auto colProc = triplet.row() + A.cols();
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

            if (initialized) {
                ZeeLogWarning << "Already applied an initial partitioning"
                    " instead refining the partitioning on A" << endLog;
                this->refine(A);
                return A;
            }

            // PHASE 1: split A in two, and construct extended matrix B
            this->initialize(A);

            auto B = TMatrix(2 * A.rows(), 2 * A.cols());
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
            using TTriplet = Zee::Triplet<TVal, TIdx>;

            // We only support IR on bi-partitionings
            if (A.procs() != 2) {
                ZeeLogError << "For now MG-IR only supports bipartitionings, "
                    "the matrix A has a " << A.procs() << "-way partitioning"
                    << endLog;
                return A;
            }

            // We store the current communication volume such that we can
            // see if it gets improved
            auto priorVolume = A.communicationVolume();
            auto& aImages = A.getImages();

            auto B = TMatrix(2 * A.rows(), 2 * A.cols());

            std::vector<TTriplet> coefficients;
            coefficients.reserve((int)(A.nonZeros() + A.rows() + A.cols()));

            // PHASE 1: We construct the matrix B according to the partitioning of A
            for (auto& t : *aImages[0]) {
                coefficients.push_back(TTriplet(
                            bRow(_phaseRow, t.row(), t.col(), A),
                            bCol(_phaseRow, t.row(), t.col(), A),
                            t.value()));
            }

            for (auto& t : *aImages[1]) {
                coefficients.push_back(TTriplet(
                            bRow(!_phaseRow, t.row(), t.col(), A),
                            bCol(!_phaseRow, t.row(), t.col(), A),
                            t.value()));
            }

            for (TIdx i = 0; i < A.cols() + A.rows(); ++i) {
                coefficients.push_back(TTriplet((TIdx)i, (TIdx)i, (TVal)1));
            }

            // We want B to be bipartitioned.
            // FIXME Support for lambda distribution.
            B.setDistributionScheme(Zee::Partitioning::cyclic, 1);
            B.setFromTriplets(coefficients.begin(), coefficients.end());
            B.spy("BTest");

            // PHASE 2: Apply KL
            // FIXME: naming is off here, why B member? only used in initialize
            // FIXME: want to keep B inbetween runs, since it does things
            KernighanLin<decltype(B)> kerLin(B);
            kerLin.run();

            // PHASE 3: retrieve partitioning for A
            inducePartitioning(A, B);

            if (A.communicationVolume() == priorVolume) {
                if (!_phaseRow && _rowOptimal) {
                    _locallyOptimal = true;
                } else {
                    if (_phaseRow)
                        _rowOptimal = true;
                    else
                        _rowOptimal = false;
                    _phaseRow = !_phaseRow;
                }
            }

            return A;
        }

    private:
        unique_ptr<vector<vector<atomic_wrapper<bool>>>> _bitInRow;
        bool initialized = false;
        bool _locallyOptimal = false;
        bool _rowOptimal = false;
        bool _phaseRow = true;
};
