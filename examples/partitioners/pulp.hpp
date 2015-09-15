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
            for(auto& image : images) {
                auto nz = image->nonZeros();
                if (nz == 0) {
                    cerr << "ERROR: Matrix image empty!" << endl;
                    continue;
                }

                for (auto& trip : *image) {
                    // actually want to run this simultaeneously for each procs
                    // think

                    // put in some kind of request
                    auto countLambda = [col_ = col, row_ = row]
                        (std::shared_ptr<typename TMatrix::image_type> _img)
                        -> typename TMatrix::index_type
                    {
                        typename TMatrix::index_type count = 0;
                        for (auto& a : *_img) {
                            if (a.col() == col_ || a.row() == row_)
                                count++;
                        }
                        return count;
                    };

                    auto counts = A.template compute<
                        typename TMatrix::index_type>(countLambda);

                    auto max = 0;
                    auto max_index = -1;
                    for (typename TMatrix::index_type i = 0;
                            i < counts.size(); ++i) {
                        if (counts[i] > max) {
                            max = counts[i];
                            max_index = i;
                        }
                    }

                    // give element to max_index
                    image->popElement(element_index);
                    images[max_index]->pushTriplet(trip);
                }
            }
        }

        void TMatrix& refine() override
        {
            refineWithIterations(default_number_iterations);
        }
};
