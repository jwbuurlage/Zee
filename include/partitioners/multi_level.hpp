// TODO:
// [ ] Row nets of size 1 can be removed
// [ ] Epsilon as param
// [ ] Fix initial sizes of A, B
// [ ] Buckets A / B, O(1) best vertex selection
// [ ] Make vertex gain updates more efficient
// [ ] Precompute column Weights
// [ ] Generalize underyling hypergraph (row, column, fine)
// [ ] What can be parallelized here?

#pragma once

#include <zee.hpp>

namespace Zee {

#include <algorithm>

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <set>
using std::set;

#include <iostream>
using std::cout;
using std::endl;

template <class TMatrix = Zee::DSparseMatrix<double>>
class MultiLevelOneD : Zee::Partitioner<TMatrix>
{
    public:
        using TIdx = typename TMatrix::index_type;
        using TImage = typename TMatrix::image_type;


        MultiLevelOneD() {
            this->procs_ = 2;
            this->procs_in_ = 1;
        }

        virtual void initialize(TMatrix& A) override
        {
            initialized = true;
            initializeHypergraph(A);
        }

        void initializeHypergraph(TMatrix& M)
        {
            // We use the row-net model since we partition by column
            // H_0 = ( V = columns, E = rows )
            // The colNets are used in particular to update 'dirty pins'

            // The row nets and column nets are initialized
            // Can be used for coarsening as well as partitioning
            rowNets.resize(M.getRows());
            colNets.resize(M.getCols());
            for (auto& pimg : M.getImages()) {
                for (auto& triplet : *pimg) {
                    rowNets[triplet.row()].push_back(triplet.col());
                    colNets[triplet.col()].push_back(triplet.row());
                }
            }
        }

        void coarsen(TMatrix& A)
        {
            ZeeLogWarning << "Multi-level coarsening not implemented yet." << endLog;
            // coarsen graph A and store HG
            //
            // We need to recursively find columns that *match* and merge them
            //
            // is KL used for partitioning or for coarsening?
        }

        vector<bool> randomColumnDistribution(TIdx cols)
        {
            // We uniformly distribute the pins over A and B
            vector<bool> columnInA(cols, 0);

            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randbool(0, 1);

            for (TIdx col = 0; col < columnInA.size(); ++col) {
                columnInA[col] = (bool)randbool(mt);
            }

            return columnInA;
        }

        vector<TIdx> obtainRowCountA(const vector<bool>& columnInA)
        {
            vector<TIdx> rowCountA(rowNets.size());
            auto r_idx = 0;
            for (auto& rowNet : rowNets) {
                for (auto& pin : rowNet) {
                    if (columnInA[pin] == true) {
                        rowCountA[r_idx]++;
                    }
                }
                r_idx++;
            }
            return rowCountA;
        }

        TIdx gainForPinInRow(bool columnInA, TIdx count, TIdx size)
        {
            if (size == 1) {
                return 0;
            } else if ((columnInA && count == 1)
                    || (!columnInA && count == size - 1)) {
                return 1;
            } else if (count == 0 || count == size) {
                return -1;
            }
            return 0;
        }

        /** For MG refine is actually a repartitioning into A^c + A^r
         * using information about the distribution of A */
        // We partition (V = ) C = A cup B
        // therefore we call the matrix M
        virtual TMatrix& partition(TMatrix& M) override
        {
            ZeeLogInfo << "Partitioning using FM heuristic" << endLog;

            ZeeLogVar(M.getProcs());

            // We check if everything has been initialized properly
            if (!initialized) {
                ZeeLogError << "Trying to partition with uninitialized partitioner." << endLog;
                return M;
            }

            // We bipartition
            TIdx p = 2;

            // The allowed load imbalance FIXME as param
            double epsilon = 0.03;

            TIdx allowedSize = (TIdx)(0.5 * M.nonZeros() * (1 + epsilon));

            // First start with a random distribution
            // This vector says for each column if it is in "A"
            vector<bool> columnInA = randomColumnDistribution(M.getCols());

            // Weights are number of elements in each column of the matrix
            // TODO: precompute this in matrix?
            vector<TIdx> columnWeights(M.getCols(), 0);
            for (TIdx j = 0; j < M.getCols(); ++j) {
                columnWeights[j] = M.getColumnWeight(j);
            }

            // We need to make sure the distribution satisfies the imbalance
            // otherwise we flip the first `x` until it does
            TIdx counts[2] = { 0, 0 };
            for (TIdx j = 0; j < M.getCols(); ++j) {
                if (columnInA[j]) {
                    counts[0] += columnWeights[j];
                } else {
                    counts[1] += columnWeights[j];
                }
            }

            if (counts[0] > allowedSize) {
                // FIXME prefix
            } else if (counts[1] > allowedSize) {
                // FIXME prefix
            }

            // Store the distribution for each net
            // rowCountA[row] yields the #elements in A in row
            vector<TIdx> rowCountA = obtainRowCountA(columnInA);

            // `max_size` is the maximum degree of a column
            // This is equal to the maximum size of a colNet
            auto max_size = std::accumulate(colNets.begin(),
                    colNets.end(), (TIdx)0,
                    [] (TIdx rhs, const vector<TIdx>& elem) {
                        return (rhs > elem.size()) ? rhs : elem.size();
                    });

            // Next we make the bucketed lists.
            // We need a vector of size 2*m + 1 containing lists
            vector<list<TIdx>> buckets(2 * max_size + 1);
            vector<TIdx> vertexGains(M.getCols());

            // Next we loop over all vertices, obtaining vertex gain
            // if it is the *only* element of {A, B} in a net, a gain of 1
            // if the net is and stays mixed then a gain of 0
            // if the net is pure then flipping it is a gain of -1
            for (TIdx row = 0; row < M.getRows(); ++row) {
                for (auto& pin : rowNets[row]) {
                    vertexGains[pin] += gainForPinInRow(columnInA[pin],
                            rowCountA[row],
                            rowNets[row].size());
                }
            }

            // We initialize the buckets, and store the list elements for each
            // column in M
            vector<typename std::list<TIdx>::iterator> listElements(M.getCols());
            for (TIdx j = 0; j < M.getCols(); ++j) {
                buckets[vertexGains[j] + max_size].push_back(j);
                listElements[j] = --buckets[vertexGains[j] + max_size].end();
            }

            // Store the pins that become 'dirty' after a flip
            set<std::pair<TIdx, TIdx>> verticesToUpdate;

            // We need to choose A or B at random if they tie for an elements
            // with largest gain
            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randbool(0, 1);

            // We copy the bool vector to store the best split
            // FIXME: instead we could just store the swaps and
            // reverse engineer
            vector<bool> bestSplit = columnInA;

            // We maintain the current and maximum net gain
            TIdx netGain = 0;
            TIdx maxNetGain = 0;

            // Choose the base cell (e.g. tail of highest non-trivial list)
            for (auto iter = 0; iter < 10000; ++iter) {
                signed long long baseCell = -1;
                // FIXME can we do this in constant time? the loop is possibly
                // unnecessary.. maybe linked list as well
                for (TIdx bucket = buckets.size() - 1; bucket >= 0; --bucket) {
                    for (auto cell : buckets[bucket]) {
                        if (counts[columnInA[cell] ? 1 : 0] >= allowedSize)
                            continue;
                        baseCell = cell;
                        netGain += bucket - max_size;
                        break;
                    }

                    if (baseCell != -1)
                        break;
                }

                // We found the baseCell
                // At the start of each iteration we need to clear the dirty vertices
                // i.e. the pins for which the vertexGains is not correct anymore
                verticesToUpdate.clear();

                if (baseCell == -1) {
                    ZeeLogError << "No viable base cell found." << endLog;
                    break;
                }

                for (auto row : colNets[baseCell]) {
                    for (auto pin : rowNets[row]) {
                        verticesToUpdate.insert(std::make_pair(pin, vertexGains[pin]));
                    }
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
                        vertexGains[pin] -= gainForPinInRow(columnInA[pin],
                                rowCountA[row],
                                rowNets[row].size());
                    }

                    rowCountA[row] += columnInA[baseCell] ? -1 : 1;

                    if (rowNets[row].size() == 1)
                        continue;

                    //// FIXME can probably be made more efficient, similar to nullifying step
                    //// case (1)
                    //if ((columnInA[baseCell] && rowCountA[row] == 1)
                    //        || (!columnInA[baseCell] && rowCountA[row] == rowNets[row].size() - 1)) {
                    //    // from positive gain to neutral gain
                    //    vertexGains[baseCell] -= 1;
                    //}
                    //// case (2)
                    //else if (rowCountA[row] == rowNets[row].size() || rowCountA[row] == 0) {
                    //    // from negative gain, to neutral gain
                    //    for (auto pin : rowNets[row]) {
                    //        vertexGains[pin] += 1;
                    //    }
                    //}
                }

                // We flip the base cell
                columnInA[baseCell] = !columnInA[baseCell];
                counts[columnInA[baseCell] ? 0 : 1] += columnWeights[baseCell];
                counts[columnInA[baseCell] ? 1 : 0] -= columnWeights[baseCell];

                // We fix the gains in the affected rows
                for (auto row : colNets[baseCell]) {
                    // FIXME can probably be made more efficient, similar to nullifying step
                    for (auto pin : rowNets[row]) {
                        vertexGains[pin] += gainForPinInRow(columnInA[pin],
                                rowCountA[row],
                                rowNets[row].size());
                    }
                }

                for (auto pPinGain : verticesToUpdate) {
                    // remove and reinsert
                    auto originalBucket = pPinGain.second + max_size;
                    auto newBucket = vertexGains[pPinGain.first] + max_size;
                    buckets[originalBucket].erase(listElements[pPinGain.first]);
                    buckets[newBucket].push_back(pPinGain.first);
                    listElements[pPinGain.first] = --buckets[newBucket].end();
                }

                // We update max gain, and copy columnInA again if it is better
                if (netGain > maxNetGain) {
                    maxNetGain = netGain;
                    bestSplit = columnInA;
                }
            }

            // Obtain partitioning on A from Hg partitioning
            vector<unique_ptr<TImage>> new_images;
            for (TIdx i = 0; i < p; ++i)
                new_images.push_back(std::make_unique<TImage>());

            for (auto& pimg : M.getImages()) {
                for (auto& triplet : *pimg) {
                    new_images[(TIdx)(bestSplit[triplet.col()])]->pushTriplet(triplet);
                }
            }

            M.resetImages(new_images);

            ZeeLogInfo << "FM partitioned into sizes: " << counts[0] << " " << counts[1] << endLog;

            return M;
        }

    private:
        bool initialized = false;

        // The nets corresponding to the hypergraph model(s)
        vector<vector<TIdx>> colNets;
        vector<vector<TIdx>> rowNets;
};

} // namespace Zee
