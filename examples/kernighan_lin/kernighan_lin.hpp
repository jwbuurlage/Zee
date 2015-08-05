template<typename TMatrix>
class Hypergraph
{
    public:
        initialize(const TMatrix& A);

        const vector<vector<TIdx>> colNets() const {
            return _colNets;
        };

        const vector<vector<TIdx>> rowNets() const {
            return _rowNets;
        };

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
        void initialize(TMatrix& A)
        {
            _A = A;

            // if bipartitioned we take this partitioning
            // as initial for hypergraph
            if (_A.procs() == 2) {
                // columnInA from A
            } else {
                // otherwise we randomize a partitioning
                _columnInA = randomColumnDistribution();
            }

            _H.initialize(_A);

            // The allowed load imbalance FIXME as param
            _epsilon = 0.03;

            _allowedSize = (TIdx)(0.5 * _A.nonZeros() * (1 + epsilon));

            // Weights are number of elements in each column of the matrix
            // TODO: precompute this in matrix?
            _columnWeights(_A.cols(), 0);
            for (TIdx j = 0; j < _A.cols(); ++j) {
                _columnWeights[j] = _A.getColumnWeight(j);
            }

            // We need to make sure the distribution satisfies the imbalance
            // otherwise we flip the first `x` until it does
            _counts[2] = { 0, 0 };
            for (TIdx j = 0; j < _A.cols(); ++j) {
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
            _maxSize = std::accumulate(H.colNets().begin(),
                    H.colNets().end(), (TIdx)0,
                    [] (TIdx rhs, const vector<TIdx>& elem) {
                        return (rhs > elem.size()) ? rhs : elem.size();
                    });

            // Next we make the bucketed lists.
            // We need a vector of size 2*m + 1 containing lists
            _buckets(2 * _maxSize + 1);
            _vertexGains(_A.cols());

            // Next we loop over all vertices, obtaining vertex gain
            // if it is the *only* element of {A, B} in a net, a gain of 1
            // if the net is and stays mixed then a gain of 0
            // if the net is pure then flipping it is a gain of -1
            for (TIdx row = 0; row < _A.rows(); ++row) {
                for (auto& pin : rowNets[row]) {
                    _vertexGains[pin] += gainForPinInRow(_columnInA[pin],
                            _rowCountA[row],
                            _H.rowNets()[row].size());
                }
            }

            // We initialize the buckets, and store the list elements for each
            // column in M
            _listElements(_A.cols());
            for (TIdx j = 0; j < _A.cols(); ++j) {
                _buckets[vertexGains[j] + _maxSize].push_back(j);
                _listElements[j] = --_buckets[vertexGains[j] + _maxSize].end();
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
            vector<TIdx> rowCountA(H.rowNets().size());
            auto r_idx = 0;
            for (auto& rowNet : H.rowNets()) {
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
        void step()
        {
            // We check if everything has been initialized properly
            if (!_initialized) {
                ZeeLogError << "Trying to partition with uninitialized partitioner." << endLog;
                return;
            }

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

            // Choose the base cell (e.g. tail of highest non-trivial list)
            for (auto iter = 0; iter < H.colNets().size; ++iter) {
                TIdx baseCell = -1;
                // FIXME can we do this in constant time? the loop is possibly
                // unnecessary.. maybe linked list as well
                for (TIdx bucket = buckets.size() - 1; bucket >= 0; --bucket) {
                    for (auto cell : buckets[bucket]) {
                        if (_counts[columnInA[cell] ? 1 : 0] >= allowedSize)
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
                _columnInA[baseCell] = !_columnInA[baseCell];
                _counts[columnInA[baseCell] ? 0 : 1] += columnWeights[baseCell];
                _counts[columnInA[baseCell] ? 1 : 0] -= columnWeights[baseCell];

                // We fix the gains in the affected rows
                for (auto row : colNets[baseCell]) {
                    // FIXME can probably be made more efficient, similar to nullifying step
                    for (auto pin : rowNets[row]) {
                        _vertexGains[pin] += gainForPinInRow(_columnInA[pin],
                                _rowCountA[row],
                                _H.rowNets()[row].size());
                    }
                }

                for (auto pPinGain : verticesToUpdate) {
                    // remove and reinsert
                    auto originalBucket = pPinGain.second + _maxSize;
                    auto newBucket = vertexGains[pPinGain.first] + _maxSize;
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

            ZeeLogInfo << "FM partitioned into sizes: " << _counts[0] << " " << _counts[1] << endLog;
        };

    private:
        // We need to preprocess the matrix before we can iterate
        bool _initialized;

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
        vector<TIdx> _columnInA;
        vector<TIdx> _columnWeights;
        vector<TIdx> _rowCountA;
        TIdx _counts[2];
        vector<list<TIdx>> _buckets;
        vector<TIdx> _vertexGains;
        vector<typename std::list<TIdx>::iterator> _listElements;
};
