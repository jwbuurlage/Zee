// we ignore distributed memory for now
// TODO:
// [ ] Fix vertex weights for imbalance

#include <zee.hpp>

#include "../hypergraph/hypergraph.hpp"

#include <random>
#include <vector>
#include <memory>

enum class HGModel { fine_grain = 1, row_net = 2, column_net = 3 };

template <typename T>
size_t argmax(std::vector<T>& list) {
    unsigned int ans = 0;
    T min = std::numeric_limits<T>::min();
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
    constexpr const static double default_load_imbalance = 0.1;

  public:
    PulpPartitioner(TMatrix& A, TIdx procs,
                    double epsilon = default_load_imbalance)
        : Zee::IterativePartitioner<TMatrix>(), A_(A) {
        this->setProcs(procs);
        // FIXME;
        minPartSize_ =
            (TIdx)(((A.nonZeros() - 1) / procs) + 1) * (1.0 - epsilon);
        ZeeLogVar(minPartSize_);
        maxPartSize_ =
            (TIdx)(((A.nonZeros() - 1) / procs) + 1) * (1.0 + epsilon);
        ZeeLogVar(maxPartSize_);
    }

    void initialize(TMatrix& A) override { initialize(A, HGModel::fine_grain); }

    void initialize(TMatrix& A, HGModel model) {
        if (initialized_)
            return;

        using TImage = typename TMatrix::image_type;

        initialized_ = true;

        if (A.getProcs() == 1) {
            std::vector<std::unique_ptr<TImage>> aNewImages;
            for (TIdx i = 0; i < this->procs_; ++i) {
                aNewImages.push_back(std::make_unique<TImage>());
            }

            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randproc(0, this->procs_ - 1);

            if (model == HGModel::fine_grain) {
                // randomly distribute A into k images
                for (auto& image : A.getImages()) {
                    for (auto& trip : *image) {
                        aNewImages[randproc(mt)]->pushTriplet(trip);
                    }
                }
            } else if (model == HGModel::row_net) {
                std::vector<TIdx> columnTargets(A.getCols());
                for (auto& target : columnTargets) {
                    target = randproc(mt);
                }

                // randomly distribute each column of A into k images
                for (auto& image : A.getImages()) {
                    for (auto& trip : *image) {
                        aNewImages[columnTargets[trip.col()]]->pushTriplet(
                            trip);
                    }
                }

            } else if (model == HGModel::column_net) {
                std::vector<TIdx> rowTargets(A.getRows());
                for (auto& target : rowTargets) {
                    target = randproc(mt);
                }

                // randomly distribute each column of A into k images
                for (auto& image : A.getImages()) {
                    for (auto& trip : *image) {
                        aNewImages[rowTargets[trip.row()]]->pushTriplet(trip);
                    }
                }
            }

            A.resetImages(aNewImages);
        } else if (A.getProcs() != this->procs_) {
            ZeeLogError << "(PuLP): Matrix is initialized with wrong number of "
                           "processors."
                        << endLog;
            return;
        }

        // FIXME should check if for row net model / column net model
        // initial matrix is distributed properly if given by user
        // make property of sparse matrix itself?

        // initialize hypergraph model?
        // yes. hypergraph model needs to be aware of Matrix
        // act like its model, but not necessarily store explicitely
        // operations we require:
        // - reassigning vertices (cached?)
        // - looping over nets
        // - apply changes to matrix

        if (model == HGModel::fine_grain)
            hyperGraph_ = std::make_unique<FineGrainHG<TIdx, TMatrix>>(A);
        else if (model == HGModel::row_net)
            hyperGraph_ = std::make_unique<RowNetHG<TIdx, TMatrix>>(A);
        else if (model == HGModel::column_net)
            hyperGraph_ = std::make_unique<ColumnNetHG<TIdx, TMatrix>>(A);

        // initalize counts
        partSize.resize(A.getProcs());
        TIdx i = 0;
        for (auto& image : A.getImages()) {
            partSize[i++] = image->nonZeros();
        }

        ZeeLogVar(hyperGraph_->getPartSize(0));
        ZeeLogVar(hyperGraph_->getPartSize(1));
    }

    virtual TMatrix& initialPartitioning(TIdx iters) {
        // hyperGraph_->sortNetsBySize();

        TIdx size = hyperGraph_->getMaximumNetSize() / 2;

        TIdx maxSize = 2;
        while (maxSize < size) {
            for (TIdx i = 0; i < iters; ++i) {
                this->refineWithIterations(A_.nonZeros(), maxSize);
            }
            ZeeLogVar(hyperGraph_->partitioningLV());
            ZeeLogVar(maxSize);
            maxSize *= 2;
        }
        return A_;
    }

    virtual TMatrix& refineWithIterations(TIdx iters, TIdx maximumNetSize = 0) {

        // initial partitioning
        TIdx r = 1;
        TIdx i = 0;
        while (r > 0 && i < iters) {
            r = 0;
            for (TIdx v = 0; v < hyperGraph_->getVertexCount(); ++v) {
                auto w = [this](TIdx peers, TIdx size) {
                    static const auto control = 0.99;
                    double z = (peers - 1.0) / (size - 1.0);
                    double x = control * (1 + z) / (1 - z);
                    return log(x);
                };

                // get neighbour counts
                auto N = hyperGraph_->partQuality(v, w, maximumNetSize);
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
