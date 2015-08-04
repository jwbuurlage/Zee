// maybe use this
template<typename TMatrix>
class Hypergraph
{
    public:
        initialize(const TMatrix& A);

    private:
        vector<vector<TIdx>> colNets;
        vector<vector<TIdx>> rowNets;
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
            // what do we need
            // 1. "bit in row" (bit in A)
            // 2. rowNets
            // 3. colNets
            //
            // if A is bipartitioned, initialize col distribution
            // from that.
            // otherwise we randomize
            if (A.procs() == 2) {

            } else {
                columnInA = randomColumnDistribution();
            }

            H.initialize(A);
        };

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

        // perform one iteration of KL
        void step();

    private:
        // bi-partitioned matrix
        // used to initialize
        vector<TIdx> columnInA;
        Hypergraph<TMatrix> H;
};
