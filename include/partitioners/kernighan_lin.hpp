#pragma once

#include <random>
#include <algorithm>
#include <vector>
#include <list>
#include <set>
#include <memory>

namespace Zee {
using std::vector;
using std::list;
using std::set;
using std::unique_ptr;

template <typename TMatrix>
class Hypergraph {
  public:
    using TIdx = typename TMatrix::index_type;

    void initialize(const TMatrix& A) {
        _rowNets.resize(A.getRows());
        _colNets.resize(A.getCols());
        for (auto& imgPtr : A.getImages()) {
            for (auto& triplet : *imgPtr) {
                _rowNets[triplet.row()].push_back(triplet.col());
                _colNets[triplet.col()].push_back(triplet.row());
            }
        }
    }

    const vector<vector<TIdx>>& colNets() const { return _colNets; }

    const vector<vector<TIdx>>& rowNets() const { return _rowNets; }

  private:
    vector<vector<TIdx>> _colNets;
    vector<vector<TIdx>> _rowNets;
};

template <typename TMatrix>
class KernighanLin {
  public:
    using TIdx = typename TMatrix::index_type;
    using TImage = typename TMatrix::image_type;

    // initialize the necessary constructs
    KernighanLin(TMatrix& A) : A_(A) {
        // if bipartitioned we take this partitioning
        // as initial for hypergraph
        if (A_.getProcs() == 2) {
            columnInA_.resize(A_.getCols());
            std::fill(columnInA_.begin(), columnInA_.end(), false);

            // columnInA from A
            // 1: get colset from A_ images
            auto colSetOne = A_.getImages()[0]->getColSet();

            // 2: set columnInA from these
            for (auto pColCount : colSetOne) {
                columnInA_[pColCount.first] = true;
            }
        } else {
            // otherwise we randomize a partitioning
            columnInA_ = randomColumnDistribution(A_.getCols());
        }

        H_.initialize(A_);

        auto& colNets = H_.colNets();
        auto& rowNets = H_.rowNets();

        // The allowed load imbalance FIXME as param
        epsilon_ = 0.03;

        allowedSize_ = (TIdx)(0.5 * A_.nonZeros() * (1 + epsilon_));

        // Weights are number of elements in each column of the matrix
        columnWeights_.resize(A_.getCols());
        for (TIdx j = 0; j < A_.getCols(); ++j) {
            columnWeights_[j] = A_.getColumnWeight(j);
        }

        // We need to make sure the distribution satisfies the imbalance
        // otherwise we flip the first `x` until it does
        counts_[0] = 0;
        counts_[1] = 0;

        for (TIdx j = 0; j < A_.getCols(); ++j) {
            if (columnInA_[j]) {
                counts_[0] += columnWeights_[j];
            } else {
                counts_[1] += columnWeights_[j];
            }
        }

        if (counts_[0] > allowedSize_) {
            // FIXME prefix
        } else if (counts_[1] > allowedSize_) {
            // FIXME prefix
        }

        // Store the distribution for each net
        // rowCountA[row] yields the #elements in A in row
        rowCountA_ = obtainRowCountA(columnInA_);

        // `max_size` is the maximum degree of a column
        // This is equal to the maximum size of a colNet
        maxSize_ =
            std::accumulate(colNets.begin(), colNets.end(), (TIdx)0,
                            [](TIdx rhs, const vector<TIdx>& elem) {
                                return (rhs > elem.size()) ? rhs : elem.size();
                            });

        // Next we make the bucketed lists.
        // We need a vector of size 2*m + 1 containing lists
        buckets_.clear();
        for (TIdx i = 0; i < 2 * maxSize_ + 1; ++i) {
            buckets_.push_back(list<TIdx>());
        }

        vertexGains_.resize(A_.getCols());
        std::fill(vertexGains_.begin(), vertexGains_.end(), 0);

        // Next we loop over all vertices, obtaining vertex gain
        // if it is the *only* element of {A, B} in a net, a gain of 1
        // if the net is and stays mixed then a gain of 0
        // if the net is pure then flipping it is a gain of -1
        for (TIdx row = 0; row < A_.getRows(); ++row) {
            for (const auto& pin : rowNets[row]) {
                vertexGains_[pin] += gainForPinInRow(
                    columnInA_[pin], rowCountA_[row], H_.rowNets()[row].size());
            }
        }

        // We initialize the buckets, and store the list elements for each
        // column in M
        listElements_.clear();
        listElements_.resize(A_.getCols());
        for (TIdx j = 0; j < A_.getCols(); ++j) {
            buckets_[vertexGains_[j] + maxSize_].push_back(j);
            listElements_[j] = --buckets_[vertexGains_[j] + maxSize_].end();
        }
    };

    vector<bool> randomColumnDistribution(TIdx cols) {
        // We uniformly distribute the pins over A and B
        vector<bool> randomColumnInA(cols, 0);

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<TIdx> randbool(0, 1);

        for (TIdx col = 0; col < randomColumnInA.size(); ++col) {
            randomColumnInA[col] = (bool)randbool(mt);
        }

        return randomColumnInA;
    }

    vector<TIdx> obtainRowCountA(const vector<bool>& columnInA) {
        vector<TIdx> rowCountA(H_.rowNets().size());
        auto r_idx = 0;
        for (auto& rowNet : H_.rowNets()) {
            for (auto& pin : rowNet) {
                if (columnInA[pin] == true) {
                    rowCountA[r_idx]++;
                }
            }
            r_idx++;
        }
        return rowCountA;
    }

    TIdx gainForPinInRow(bool columnInA, TIdx count, TIdx size) {
        if (size == 1) {
            return 0;
        } else if ((columnInA && count == 1) ||
                   (!columnInA && count == size - 1)) {
            return 1;
        } else if (count == 0 || count == size) {
            return -1;
        }
        return 0;
    }

    // perform one iteration of KL
    void run() {
        // We check if everything has been initialized properly
        // we bipartition
        TIdx p = 2;

        // FIXME: We reset the partitioning
        // to the current bipartitioning of A_

        // We need to choose A or B at random if they tie for an elements
        // with largest gain
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<TIdx> randbool(0, 1);

        // We copy the bool vector to store the best split
        // FIXME: instead we could just store the swaps and
        // reverse engineer
        vector<bool> bestSplit = columnInA_;

        // We maintain the current and maximum net gain
        TIdx netGain = 0;
        TIdx maxNetGain = 0;

        // Store the pins that become 'dirty' after a flip
        set<std::pair<TIdx, TIdx>> verticesToUpdate;

        auto& colNets = H_.colNets();
        auto& rowNets = H_.rowNets();

        // Choose the base cell (e.g. tail of highest non-trivial list)
        for (TIdx iter = 0; iter < H_.colNets().size(); ++iter) {
            // FIXME: we want to change this to support BigNums etc.
            signed long baseCell = -1;

            // FIXME can we do this in constant time? the loop is possibly
            // unnecessary.. maybe linked list as well
            for (TIdx bucket = buckets_.size() - 1; bucket >= 0; --bucket) {
                for (auto cell : buckets_[bucket]) {
                    if (counts_[columnInA_[cell] ? 1 : 0] + 1 > allowedSize_)
                        continue;
                    baseCell = cell;
                    netGain += bucket - maxSize_;
                    break;
                }

                if (baseCell != -1)
                    break;
            }

            // We found the baseCell
            // At the start of each iteration we need to clear the dirty
            // vertices
            // i.e. the pins for which the vertexGains is not correct anymore
            verticesToUpdate.clear();

            if (baseCell == -1) {
                ZeeLogError << "No viable base cell found." << endLog;
                break;
            }

            // If rowCountA changes 'type', then we need to update pins
            // The types are:
            // (1) 'almost' A or B (1 or size - 1)
            // (2) 'completely' A or B (0 or size)
            // (3) mixed
            // Note: if we are now at (2), we *must* have come from (1)
            // but if we got to (1) we *might* have gotten from (2)
            //
            // first we nullify the effects of the effected nets (if necessary)
            for (auto row : colNets[baseCell]) {
                for (auto pin : rowNets[row]) {
                    vertexGains_[pin] -=
                        gainForPinInRow(columnInA_[pin], rowCountA_[row],
                                        H_.rowNets()[row].size());
                }

                // FIXME: do we require this?
                //  if (H_.rowNets()[row].size() == 1)
                //      continue;

                rowCountA_[row] += columnInA_[baseCell] ? -1 : 1;
            }

            // We flip the base cell
            columnInA_[baseCell] = !columnInA_[baseCell];
            counts_[columnInA_[baseCell] ? 0 : 1] += columnWeights_[baseCell];
            counts_[columnInA_[baseCell] ? 1 : 0] -= columnWeights_[baseCell];

            // We fix the gains in the affected rows
            for (auto row : H_.colNets()[baseCell]) {
                // FIXME can probably be made more efficient, similar to
                // nullifying step
                for (auto pin : H_.rowNets()[row]) {
                    vertexGains_[pin] +=
                        gainForPinInRow(columnInA_[pin], rowCountA_[row],
                                        H_.rowNets()[row].size());
                }
            }

            for (auto row : H_.colNets()[baseCell]) {
                for (auto pin : H_.rowNets()[row]) {
                    verticesToUpdate.insert(
                        std::make_pair(pin, vertexGains_[pin]));
                }
            }

            for (auto pPinGain : verticesToUpdate) {
                // remove and reinsert
                auto originalBucket = pPinGain.second + maxSize_;
                auto newBucket = vertexGains_[pPinGain.first] + maxSize_;
                buckets_[originalBucket].erase(listElements_[pPinGain.first]);
                buckets_[newBucket].push_back(pPinGain.first);
                listElements_[pPinGain.first] = --buckets_[newBucket].end();
            }

            // We update max gain, and copy columnInA again if it is better
            if (netGain > maxNetGain) {
                maxNetGain = netGain;
                bestSplit = columnInA_;
            }
        }

        // Obtain partitioning on A from Hg partitioning
        vector<unique_ptr<TImage>> newImages;
        for (TIdx i = 0; i < p; ++i)
            newImages.push_back(std::make_unique<TImage>());

        for (auto& pimg : A_.getImages()) {
            for (auto& triplet : *pimg) {
                newImages[(TIdx)(bestSplit[triplet.col()])]->pushTriplet(
                    triplet);
            }
        }

        A_.resetImages(newImages);
    };

  private:
    // Bi-partitioned matrix
    TMatrix& A_;

    // Underyling hypergraph of A
    Hypergraph<TMatrix> H_;

    // The maximum gain of a vertex, equal to the maximum row degree
    TIdx maxSize_;

    // The 'load-imbalance'
    double epsilon_;
    TIdx allowedSize_;

    // We cache various information on the hypergraph
    vector<bool> columnInA_;
    vector<TIdx> columnWeights_;
    vector<TIdx> rowCountA_;
    TIdx counts_[2];
    vector<list<TIdx>> buckets_;
    vector<TIdx> vertexGains_;
    vector<typename std::list<TIdx>::iterator> listElements_;
};

} // namespace Zee
