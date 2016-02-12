#pragma once

#include <zee.hpp>

#include "../hypergraph/hypergraph.hpp"

#include <random>
#include <vector>
#include <memory>
#include <limits>

namespace Zee {

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
    using TImage = typename TMatrix::image_type;

    constexpr const static TIdx default_number_iterations = 100;
    constexpr const static double default_load_imbalance = 0.1;

  public:
    PulpPartitioner(TMatrix& A, TIdx procs,
                    double epsilon = default_load_imbalance)
        : Zee::IterativePartitioner<TMatrix>(), A_(A) {
        this->setProcs(procs);
        minPartSize_ =
            (TIdx)(((A.nonZeros() - 1) / procs) + 1) * (1.0 - epsilon);
        maxPartSize_ =
            (TIdx)(((A.nonZeros() - 1) / procs) + 1) * (1.0 + epsilon);
        epsilon_ = epsilon;
    }

    void initialize(TMatrix& A) override { initialize(A, HGModel::fine_grain); }

    void randomReset(TMatrix& A, HGModel model) {
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
                    aNewImages[columnTargets[trip.col()]]->pushTriplet(trip);
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


        constructHypergraph(A, model);


        if (A.loadImbalance() > (1.0 + epsilon_)) {
            if (this->procs_ > 2) {
                ZeeLogError << "Invalid initial size for non-bipartitioning. "
                               "No support for fixing the initialization for "
                               "small matrices."
                            << endLog;
                exit(-1);
            }
            TIdx largePart = 0;
            if (hyperGraph_->getPartSize(1) > hyperGraph_->getPartSize(0)) {
                    largePart = 1;
            }

            auto weights = hyperGraph_->getWeights();
            std::vector<TIdx> X(hyperGraph_->getVertexCount());
            std::iota(X.begin(), X.end(), 0);
            std::sort(X.begin(), X.end(), [weights](TIdx lhs, TIdx rhs) {
                return weights[lhs] < weights[rhs];
            });

            for (TIdx i = 0; i < hyperGraph_->getVertexCount(); ++i) {
                auto v = X[i];
                if (hyperGraph_->getPart(v) == largePart)
                    hyperGraph_->reassign(v, 1 - largePart);
                if (A.loadImbalance() < (1.0 + epsilon_))
                    break;
            }
            if (A.loadImbalance() > (1.0 + epsilon_)) {
                ZeeLogError << "Cannot find valid initial distribution. eps = "
                            << A.loadImbalance() << endLog;
                exit(-1);
            }
        }
    }

    void constructHypergraph(TMatrix& A, HGModel model) {
        if (model == HGModel::fine_grain)
            hyperGraph_ = std::make_unique<FineGrainHG<TIdx, TMatrix>>(A);
        else if (model == HGModel::row_net)
            hyperGraph_ = std::make_unique<RowNetHG<TIdx, TMatrix>>(A);
        else if (model == HGModel::column_net)
            hyperGraph_ = std::make_unique<ColumnNetHG<TIdx, TMatrix>>(A);
    }

    void initialize(TMatrix& A, HGModel model) {
        if (initialized_)
            return;
        initialized_ = true;

        // FIXME should check if for row net model / column net model
        // initial matrix is distributed properly if given by user
        // make property of sparse matrix itself?

        if (A.getProcs() == 1) {
            randomReset(A, model);
        } else {
            constructHypergraph(A, model);
        }
    }

    virtual TMatrix& initialPartitioning(TIdx iters) {
        hyperGraph_->sortNetsBySize();

        TIdx maximumNetCount = hyperGraph_->getNetCount() / 2;

        TIdx netsToConsider = 1;
        while (netsToConsider < maximumNetCount) {
            for (TIdx i = 0; i < iters; ++i) {
                this->refineWithIterations(A_.nonZeros(), 0, netsToConsider);
            }
            netsToConsider *= 2;
        }
        return A_;
    }

    virtual TMatrix& refineWithIterations(TIdx iters, TIdx maximumNetSize = 0,
                                          TIdx netsToConsider = 0, bool randomize = false) {

        // initial partitioning
        TIdx r = 1;
        TIdx i = 0;
        while (r > 0 && i < iters) {
            r = 0;

            std::random_device rd;
            std::mt19937 mt(rd());

            std::vector<TIdx> X;
            if (randomize) {
                X.resize(hyperGraph_->getVertexCount());
                std::iota(X.begin(), X.end(), 0);
                std::shuffle(X.begin(), X.end(), mt);
            }

            for (TIdx j = 0; j < hyperGraph_->getVertexCount(); ++j) {
                auto v = j;
                if (randomize)
                    v = X[j];

                auto w = [](double peers, double size) {
                    static const auto control = 0.99;
                    double z = control * (((2.0 * peers) / size) - 1.0);
                    double x = (1 + z) / (1 - z);
                    return log(x);
                };

                // get neighbour counts
                auto N = hyperGraph_->partQuality(v, w, maximumNetSize,
                                                  netsToConsider);
                auto p = argmax(N);

                if (p != hyperGraph_->getPart(v) &&
                    hyperGraph_->getPartSize(p) +
                            hyperGraph_->getWeights()[v] <=
                        maxPartSize_) {
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

    double epsilon_;
    TIdx minPartSize_;
    TIdx maxPartSize_;

    std::unique_ptr<DHypergraph<TIdx>> hyperGraph_;

    // neighbour histogram (N_v)
    std::vector<std::vector<TIdx>> neighbourCount;
};

} // namespace Zee
