#include <zee.hpp>
#include <random>

#include <memory>
using std::shared_ptr;


template <class TMatrix = Zee::DSparseMatrix<>>
class PulpPartitioner : Zee::IterativePartitioner<TMatrix>
{
    using TIdx = TMatrix::index_type;
    //using TVal = TMatrix::value_type;
    constexpr const static TIdx default_number_iterations = 100;

    public:
        PulpPartitioner(TMatrix& A, TIdx procs) :
            Zee::IterativePartitioner(A)
        {
            setProcs(procs);
            initialize();
        } 

        void initialize() override {
            if (initialized_)
                return;

            // initalize structures
        }

        virtual TMatrix& refineWithIterations(TIdx iters) {
            auto& images = A.getMutableImages();
            auto p = A.getProcs();

            void distCs = [] (
                std::shared_ptr<typename TMatrix::image_type> _img,
                TIdx s)
            {
                for (auto& pRowCount : img->getRowSet()) {
                    auto row = pRowCount.first;
                    auto count = pRowCount.second;
                    
                    // TODO Fixme syntax how
                    Zee::message::put(s, row % p);
                    Zee::message::put(row, row % p);
                    Zee::message::put(count, row % p);
                }

                for (auto& pColCount : img->getColSet()) {
                    auto col = pColCount.first;
                    auto count = pColCount.second;
                   
                    // TODO Fixme
                    Zee::message::put(s, col % p);
                    Zee::message::put(row, col % p);
                    Zee::message::put(count, col % p);
                }
            };

            void redistCs = [] (
                std::shared_ptr<typename TMatrix::image_type> _img,
                TIdx s)
            {

            };
        }

        void TMatrix& refine() override
        {
            refineWithIterations(default_number_iterations);
        }
};
