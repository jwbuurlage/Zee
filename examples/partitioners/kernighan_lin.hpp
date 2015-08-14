#include <random>
#include <algorithm>

#include <vector>
using std::vector;

#include <list>
using std::list;

#include <set>
using std::set;

#include <memory>
using std::unique_ptr;

template<typename TMatrix>
class Hypergraph
{
    public:
        using TIdx = typename TMatrix::index_type;

        void initialize(const TMatrix& A)
        {
            _rowNets.resize(A.getRows());
            _colNets.resize(A.getCols());
            for (auto& imgPtr : A.getImages()) {
                for (auto& triplet : *imgPtr) {
                    _rowNets[triplet.row()].push_back(triplet.col());
                    _colNets[triplet.col()].push_back(triplet.row());
                }
            }
        }

        const vector<vector<TIdx>>& colNets() const {
            return _colNets;
        }

        const vector<vector<TIdx>>& rowNets() const {
            return _rowNets;
        }

    private:
        vector<vector<TIdx>> _colNets;
        vector<vector<TIdx>> _rowNets;
};

template<typename TMatrix>
class KernighanLin
{
    public:
        using TIdx = typename TMatrix::index_type;
        using TImage = typename TMatrix::image_type;

        // initialize the necessary constructs
        KernighanLin(TMatrix& A) : _A(A)
        {
            // if bipartitioned we take this partitioning
            // as initial for hypergraph
            if (_A.getProcs() == 2) {
                _columnInA.resize(_A.getCols());
                std::fill(_columnInA.begin(), _columnInA.end(), false);

                // columnInA from A
                // 1: get colset from _A images
                auto colSetOne = _A.getImages()[0]->getColSet();

                // 2: set columnInA from these
                for (auto pColCount : colSetOne) {
                    _columnInA[pColCount.first] = true;
                }
            } else {
                // otherwise we randomize a partitioning
                _columnInA = randomColumnDistribution(_A.getCols());
            }

            _H.initialize(_A);

            auto& colNets = _H.colNets();
            auto& rowNets = _H.rowNets();

            // The allowed load imbalance FIXME as param
            _epsilon = 0.03;

            _allowedSize = (TIdx)(0.5 * _A.nonZeros() * (1 + _epsilon));

            // Weights are number of elements in each column of the matrix
            _columnWeights.resize(_A.getCols());
            for (TIdx j = 0; j < _A.getCols(); ++j) {
                _columnWeights[j] = _A.getColumnWeight(j);
            }

            // We need to make sure the distribution satisfies the imbalance
            // otherwise we flip the first `x` until it does
            _counts[0] = 0;
            _counts[1] = 0;

            for (TIdx j = 0; j < _A.getCols(); ++j) {
                if (_columnInA[j]) {
                    _counts[0] += _columnWeights[j];
                } else {
                    _counts[1] += _columnWeights[j];
                }
            }

            if (_counts[0] > _allowedSize) {
                // FIXME prefix
            } else if (_counts[1] > _allowedSize) {
                // FIXME prefix
            }

            // Store the distribution for each net
            // rowCountA[row] yields the #elements in A in row
            _rowCountA = obtainRowCountA(_columnInA);

            // `max_size` is the maximum degree of a column
            // This is equal to the maximum size of a colNet
            _maxSize = std::accumulate(colNets.begin(),
                    colNets.end(), (TIdx)0,
                    [] (TIdx rhs, const vector<TIdx>& elem) {
                        return (rhs > elem.size()) ? rhs : elem.size();
                    });

            // Next we make the bucketed lists.
            // We need a vector of size 2*m + 1 containing lists
            _buckets.clear();
            for (TIdx i = 0; i < 2 * _maxSize + 1; ++i) {
                _buckets.push_back(list<TIdx>());
            }

            _vertexGains.resize(_A.getCols());
            std::fill(_vertexGains.begin(), _vertexGains.end(), 0);

            // Next we loop over all vertices, obtaining vertex gain
            // if it is the *only* element of {A, B} in a net, a gain of 1
            // if the net is and stays mixed then a gain of 0
            // if the net is pure then flipping it is a gain of -1
            for (TIdx row = 0; row < _A.getRows(); ++row) {
                for (const auto& pin : rowNets[row]) {
                    _vertexGains[pin] += gainForPinInRow(_columnInA[pin],
                            _rowCountA[row],
                            _H.rowNets()[row].size());
                }
            }

            // We initialize the buckets, and store the list elements for each
            // column in M
            _listElements.clear();
            _listElements.resize(_A.getCols());
            for (TIdx j = 0; j < _A.getCols(); ++j) {
                _buckets[_vertexGains[j] + _maxSize].push_back(j);
                _listElements[j] = --_buckets[_vertexGains[j] + _maxSize].end();
            }

        };

        vector<bool> randomColumnDistribution(TIdx cols)
        {
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

        vector<TIdx> obtainRowCountA(const vector<bool>& columnInA)
        {
            vector<TIdx> rowCountA(_H.rowNets().size());
            auto r_idx = 0;
            for (auto& rowNet : _H.rowNets()) {
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

        // perform one iteration of KL
        void run()
        {
            // We check if everything has been initialized properly
            // we bipartition
            auto p = 2;

            // FIXME: We reset the partitioning
            // to the current bipartitioning of _A

            // We need to choose A or B at random if they tie for an elements
            // with largest gain
            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randbool(0, 1);

            // We copy the bool vector to store the best split
            // FIXME: instead we could just store the swaps and
            // reverse engineer
            vector<bool> bestSplit = _columnInA;

            // We maintain the current and maximum net gain
            TIdx netGain = 0;
            TIdx maxNetGain = 0;

            // Store the pins that become 'dirty' after a flip
            set<std::pair<TIdx, TIdx>> verticesToUpdate;

            auto& colNets = _H.colNets();
            auto& rowNets = _H.rowNets();

            // Choose the base cell (e.g. tail of highest non-trivial list)
            for (auto iter = 0; iter < _H.colNets().size(); ++iter) {
                TIdx baseCell = -1;

                // FIXME can we do this in constant time? the loop is possibly
                // unnecessary.. maybe linked list as well
                for (TIdx bucket = _buckets.size() - 1; bucket >= 0; --bucket) {
                    for (auto cell : _buckets[bucket]) {
                        if (_counts[_columnInA[cell] ? 1 : 0] + 1 > _allowedSize)
                            continue;
                        baseCell = cell;
                        netGain += bucket - _maxSize;
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
                        _vertexGains[pin] -= gainForPinInRow(_columnInA[pin],
                                _rowCountA[row],
                                _H.rowNets()[row].size());
                    }

                    // FIXME: do we require this?
                    //  if (_H.rowNets()[row].size() == 1)
                    //      continue;

                    _rowCountA[row] += _columnInA[baseCell] ? -1 : 1;
                }

                // We flip the base cell
                _columnInA[baseCell] = !_columnInA[baseCell];
                _counts[_columnInA[baseCell] ? 0 : 1] += _columnWeights[baseCell];
                _counts[_columnInA[baseCell] ? 1 : 0] -= _columnWeights[baseCell];

                // We fix the gains in the affected rows
                for (auto row : _H.colNets()[baseCell]) {
                    // FIXME can probably be made more efficient, similar to nullifying step
                    for (auto pin : _H.rowNets()[row]) {
                        _vertexGains[pin] += gainForPinInRow(_columnInA[pin],
                                _rowCountA[row],
                                _H.rowNets()[row].size());
                    }
                }

                for (auto row : _H.colNets()[baseCell]) {
                    for (auto pin : _H.rowNets()[row]) {
                        verticesToUpdate.insert(std::make_pair(pin, _vertexGains[pin]));
                    }
                }

                for (auto pPinGain : verticesToUpdate) {
                    // remove and reinsert
                    auto originalBucket = pPinGain.second + _maxSize;
                    auto newBucket = _vertexGains[pPinGain.first] + _maxSize;
                    _buckets[originalBucket].erase(_listElements[pPinGain.first]);
                    _buckets[newBucket].push_back(pPinGain.first);
                    _listElements[pPinGain.first] = --_buckets[newBucket].end();
                }

                // We update max gain, and copy columnInA again if it is better
                if (netGain > maxNetGain) {
                    maxNetGain = netGain;
                    bestSplit = _columnInA;
                }
            }

            // Obtain partitioning on A from Hg partitioning
            vector<unique_ptr<TImage>> newImages;
            for (TIdx i = 0; i < p; ++i)
                newImages.push_back(std::make_unique<TImage>());

            for (auto& pimg : _A.getImages()) {
                for (auto& triplet : *pimg) {
                    newImages[(TIdx)(bestSplit[triplet.col()])]->pushTriplet(triplet);
                }
            }

            _A.resetImages(newImages);
        };

    private:
        // Bi-partitioned matrix
        TMatrix& _A;

        // Underyling hypergraph of A
        Hypergraph<TMatrix> _H;

        // The maximum gain of a vertex, equal to the maximum row degree
        TIdx _maxSize;

        // The 'load-imbalance'
        double _epsilon;
        TIdx _allowedSize;

        // We cache various information on the hypergraph
        vector<bool> _columnInA;
        vector<TIdx> _columnWeights;
        vector<TIdx> _rowCountA;
        TIdx _counts[2];
        vector<list<TIdx>> _buckets;
        vector<TIdx> _vertexGains;
        vector<typename std::list<TIdx>::iterator> _listElements;
};
