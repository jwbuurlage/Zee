#include <zee.hpp>

template <class TMatrix = Zee::DSparseMatrix<double>>
class MultiLevelOneD : Zee::Partitioner<TMatrix>
{
    public:
        MultiLevelOneD() {
            this->_procs = 2;
            this->_procs_in = 1;
        }

        virtual void initialize(TMatrix& A) override
        {
            initialized = true;

            // do we coarsen graph here?
            // lets just explicitely use KL, maybe look up paper that specializes
            // it to HGs
            //
            // We use the row-net model since we partition by column
            // H_0 = ( V = columns, E = rows )
            // We need to recursively find columns that match and merge them
            //
            // is KL used for partitioning or for coarsening?
            //
            // ignore coarsening
        }

        /** For MG refine is actually a repartitioning into A^c + A^r
         * using information about the distribution of A */
        virtual TMatrix& partition(TMatrix& A) override
        {
            // FIXME: many of these operations can be done in place

            if (!initialized) {
                Zee::logError("Trying to partition with uninitialized partitioner");
            }

            using TIdx = typename TMatrix::index_type;
            using TVal = typename TMatrix::value_type;
            using TImage = typename TMatrix::image_type;

            auto p = this->_procs;

            // perform KL on HG
            // Start with random partitioning of HG
            //
            // lets obtain the row-nets centrally for now
            vector<vector<TIdx>> rowNets(A.rows());

            vector<TIdx> columnInA(A.nonZeros());
            
            // obtain partitioning on A from HG

            return A;
        }

        void coarsen(TMatrix& A)
        {
            // coarsen graph A and store HG
        }

    private:
        bool initialized = false;
};
