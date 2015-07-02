#pragma once

#include <zee.hpp>

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
        virtual TMatrix& partition(TMatrix& M) override
        {
            ZeeInfoLog << "Partitioning using FM heuristic" << endLog;

            if (!initialized) {
                ZeeInfoLog << "Trying to partition with uninitialized partitioner" << endLog;
            }

            using TIdx = typename TMatrix::index_type;
            using TVal = typename TMatrix::value_type;
            using TImage = typename TMatrix::image_type;

            auto p = 2;

            // perform KL on HG
            // Start with random partitioning of HG
            //
            // lets obtain the row-nets centrally for now
            // and also the rows for each column (identical to col-nets)
            vector<vector<TIdx>> rowNets(M.rows());
            vector<vector<TIdx>> colNets(M.cols());
            for (auto& pimg : M.getImages()) {
                for (auto& triplet : *pimg) {
                    rowNets[triplet.row()].push_back(triplet.col());
                    colNets[triplet.col()].push_back(triplet.row());
                }
            }


            // our heuristic implementation revolves around a bucketlist for the range
            // ( -m, ..., +m )
            // where m = max_i n_i
            //
            // We partition (V = ) C = A cup B
            // therefore we call the matrix M

            // first start with a random distribution
            // this vector says for each column if it is in "A"
            vector<TIdx> columnInA(M.cols(), 0);

            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randbool(0, 1);

            for (auto& col : columnInA) {
                col = randbool(mt);
            }

            // store the degree of each vertex
            vector<TIdx> columnDegrees(M.cols());

            // store the distribution for each net 
            // rND[row] = #elements in A
            vector<TIdx> rowNetDistribution(rowNets.size());
            auto r_idx = 0;
            for (auto& rowNet : rowNets) {
                for (auto& pin : rowNet) {
                    if (columnInA[pin] == 1) {
                        rowNetDistribution[r_idx]++;
                    }
                    columnDegrees[pin]++;
                }
                r_idx++;
            }

            // max_size is the maximum degree of a column
            auto max_size = *std::max_element(columnDegrees.begin(),
                    columnDegrees.end());

            // now lets make the bucketed lists.
            // We need a vector of size 2*m + 1 containing lists
            vector<list<TIdx>> buckets(2 * max_size + 1);

            // next we loop over all vertices, obtaining vertex gain
            // if it is the *only* element of {A, B} in a net, a gain of 1
            // if the net is and stays mixed then a gain of 0
            // if the net is pure then flipping it is a gain of -1
            vector<TIdx> vertexGain(M.cols());
            for (TIdx i = 0; i < M.cols(); ++i) {
                for (auto& pin : rowNets[i]) {
                    if (columnInA[pin] && rowNetDistribution[i] == 1)
                        vertexGain[pin] += 1;
                    else if (!columnInA[pin] &&
                            rowNetDistribution[i] == rowNets[i].size() - 1)
                        vertexGain[pin] += 1;
                    else if (rowNetDistribution[i] == 0 ||
                            rowNetDistribution[i] == rowNets[i].size())
                        vertexGain[pin] -= 1;
                }
            }

            // let's see if this vertex gain shit is correct
//            for (auto& rowNet : rowNets) {
//                cout << "[ ";
//                for (auto& col : rowNet) {
//                    cout << "(" << col << ", " << columnInA[col] << "), ";
//                }
//                cout << "]" << endl;
//            }
//
//            for (TIdx j = 0; j < M.cols(); ++j) {
//                cout << j << ", " << vertexGain[j] << endl;
//            }
            // seems to be..

            vector<typename std::list<TIdx>::iterator> listElements(M.cols());
            for (TIdx j = 0; j < M.cols(); ++j) {
                buckets[vertexGain[j] + max_size].push_back(j);
                listElements[j] = --buckets[vertexGain[j] + max_size].end();
            }
            
//            auto bucket_idx  = -max_size;
//            for (auto& bucket : buckets) {
//                cout << bucket_idx++ << " = [ ";
//                for (auto& col : bucket) {
//                    cout << col << ", ";
//                }
//                cout << "]" << endl;
//            }

            // choose the base cell (e.g. tail of highest non-trivial list)
            for (auto iter = 0; iter < 4; ++iter) {
                auto base_cell = -1;
                for (auto bucket = buckets.size() - 1; bucket >= 0; --bucket) {
                    if (!buckets[bucket].empty()) {
                        base_cell = buckets[bucket].back();
                        break;
                    }
                }

                // terminate if no more cells are found
                if (base_cell == -1)
                    break;

                // for each net belonging to base_cell ..
                // and update other cells' gain
                for (auto& row : colNets[base_cell]) {
                    // need to update these rows
                    // first update distribution
                    rowNetDistribution[row] +=
                        columnInA[base_cell] ? -1 : 1;
                    
                    // then for each element in row we need to update
                    // need ptrs
                    // got 'em
                }

                // we flip base_cell
                columnInA[base_cell] = !columnInA[base_cell];
            }


            // obtain partitioning on A from HG
            vector<unique_ptr<TImage>> new_images;
            for (TIdx i = 0; i < p; ++i)
                new_images.push_back(std::make_unique<TImage>());
            
            for (auto& pimg : M.getImages()) {
                for (auto& triplet : *pimg) {
                    new_images[columnInA[triplet.col()]]->pushTriplet(triplet);
                }
            }

            M.resetImages(new_images);

            return M;
        }

        void coarsen(TMatrix& A)
        {
            // coarsen graph A and store HG
        }

    private:
        bool initialized = false;
};
