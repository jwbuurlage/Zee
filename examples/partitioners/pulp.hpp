// we ignore distributed memory for now

#include <zee.hpp>

#include "../hypergraph/hypergraph.hpp"

#include <random>
#include <vector>
#include <memory>


// positive number only
template <typename T>
size_t argmax(std::vector<T>& list) {
    unsigned int ans = 0;
    T min = 0;
    for (unsigned int i = 0; i < list.size(); ++i) {
        if (list[i] > min) {
            min = list[i];
            ans = i;
        }
    }

    return ans;
}


template <class TMatrix = Zee::DSparseMatrix<>>
class PulpPartitioner : Zee::IterativePartitioner<TMatrix> {
    using TIdx = typename TMatrix::index_type;
    using TVal = typename TMatrix::value_type;

    constexpr const static TIdx default_number_iterations = 100;

  public:
    PulpPartitioner(TMatrix& A, TIdx procs)
        : Zee::IterativePartitioner<TMatrix>(), A_(A) {
        this->setProcs(procs);
        initialize(A);
        // FIXME;
        minPartSize_ = 50;
        maxPartSize_ = 300;
    }

    void initialize(TMatrix& A) override {
        if (initialized_)
            return;

        using TImage = typename TMatrix::image_type;

        initialized_ = true;

        if (A.getProcs() == 1) {
            // randomly distribute A into k images
            std::vector<std::unique_ptr<TImage>> aNewImages;
            for (TIdx i = 0; i < this->procs_; ++i) {
                aNewImages.push_back(std::make_unique<TImage>());
            }

            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randproc(0, this->procs_ - 1);

            for (auto& image : A.getImages()) {
                for (auto& trip : *image) {
                    aNewImages[randproc(mt)]->pushTriplet(trip);
                }
            }

            A.resetImages(aNewImages);
        } else if (A.getProcs() != this->procs_) {
            ZeeLogError << "(PuLP): Matrix is initialized with wrong number of "
                           "processors."
                        << endLog;
            return;
        }

        // initialize hypergraph model?
        // yes. hypergraph model needs to be aware of Matrix
        // act like its model, but not necessarily store explicitely
        // operations we require:
        // - reassigning vertices (cached?)
        // - looping over nets
        // - apply changes to matrix

        // initalize counts
        partSize.resize(A.getProcs());
        TIdx i = 0;
        for (auto& image : A.getImages()) {
            partSize[i++] = image->nonZeros();
        }

        hyperGraph_ =
            std::make_unique<FineGrainHG<TIdx, TMatrix>>(A);
    }

    virtual TMatrix& refineWithIterations(TIdx iters) {

        // initial partitioning
        TIdx r = 1;
        TIdx i = 0;
        while (r > 0 && i < iters) {
            r = 0;
            for (TIdx v = 0; v < hyperGraph_->getVertexCount(); ++v) {
                // get neighbour counts
                auto N = hyperGraph_->getNeighbourDistribution(v);
                auto p = argmax(N);

                if (p != hyperGraph_->getPart(v) &&
                    hyperGraph_->getPartSize(hyperGraph_->getPart(v)) >
                        minPartSize_) {
                    hyperGraph_->reassign(v, p);
                    r += 1;
                }

                i++;
                if (i >= iters)
                    break;
            }
        }

        return A_;
    }

    TMatrix& refine(TMatrix& A) override {
        refineWithIterations(default_number_iterations);
        return A_;
    }

  private:
    bool initialized_ = false;
    TMatrix& A_;

    TIdx minPartSize_;
    TIdx maxPartSize_;

    std::unique_ptr<DHypergraph<TIdx>> hyperGraph_;

    // counts of parts (C)
    std::vector<TIdx> partSize;

    // neighbour histogram (N_v)
    std::vector<std::vector<TIdx>> neighbourCount;
};
