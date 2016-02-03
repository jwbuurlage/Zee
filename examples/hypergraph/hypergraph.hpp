#include <zee.hpp>

#include <vector>
#include <algorithm>

template <typename TIdx = Zee::default_index_type>
class DHypergraph {
  public:
    DHypergraph(TIdx vertexCount, TIdx partCount)
        : vertexCount_(vertexCount), partCount_(partCount),
          partSize_(partCount), part_(vertexCount), weights_(vertexCount) {}

    virtual std::vector<double> partQuality(TIdx v,
                                            std::function<double(TIdx, TIdx)> w,
                                            unsigned int maximumNetSize = 0,
                                            unsigned int netsToConsider = 0) {
        if (maximumNetSize == 0) {
            maximumNetSize = getMaximumNetSize();
        }
        if (netsToConsider == 0) {
            netsToConsider = nets_.size();
        }

        std::vector<double> result(this->partCount_);

        std::vector<TIdx> members(this->partCount_);
        for (auto n : this->netsForVertex_[v]) {
            if (this->nets_[n].size() > maximumNetSize)
                continue;

            if (this->netsSizeRank_[n] >= netsToConsider)
                continue;

            for (TIdx i = 0; i < this->partCount_; ++i)
                members[i] = 0;

            for (auto u : this->nets_[n]) {
                members[this->part_[u]]++;
            }

            for (TIdx i = 0; i < this->partCount_; ++i) {
                result[i] += w(members[i], this->nets_[n].size());
            }
        }

        return result;
    }

    const std::vector<TIdx>& getWeights() const { return weights_; }

    TIdx partitioningLV() const {
        // compute (l - 1)-metric
        TIdx result = 0;
        for (TIdx net = 0; net < netCount_; ++net) {
            TIdx lambda = 0;
            for (auto cnt : netDistribution_[net])
                if (cnt > 0)
                    lambda++;
            result += lambda - 1;
        }

        return result;
    }

    std::vector<std::vector<TIdx>>& getNets() const { return nets_; }
    virtual void reassign(TIdx vertex, TIdx part) = 0;

    void sortNetsBySize() {
        this->netsSizeRank_.resize(nets_.size());

        std::vector<TIdx> netsBySize(nets_.size());
        std::iota(netsBySize.begin(), netsBySize.end(), 0);

        std::sort(netsBySize.begin(), netsBySize.end(),
                  [this](TIdx lhs,
                          TIdx rhs) {
                      return this->nets_[lhs].size() < this->nets_[rhs].size();
                  });

        // for (auto net : netsBySize) {
        //     ZeeLogVar(this->nets_[net].size());
        // }
        // ZeeAssert(0);

        for (TIdx i = 0; i < nets_.size(); ++i) {
            this->netsSizeRank_[netsBySize[i]] = i;
        }
    }

    TIdx getMaximumNetSize() {
        if (maximumNetSize_ == 0) {
            for (auto& net : nets_) {
                if (net.size() > maximumNetSize_)
                    maximumNetSize_ = net.size();
            }
        }
        return maximumNetSize_;
    }

    TIdx getVertexCount() const { return vertexCount_; }
    TIdx getNetCount() const { return netCount_; }

    TIdx getPart(TIdx v) const { return part_[v]; }
    TIdx getPartSize(TIdx p) const { return partSize_[p]; }

  protected:
    // number of vertices in this hypergraph
    TIdx vertexCount_;

    // number of nets in this hypergraph
    TIdx netCount_;

    // maximum size of a net
    TIdx maximumNetSize_ = 0;

    // nets indices by size
    std::vector<TIdx> netsSizeRank_;

    // number of parts for this hypergraph
    TIdx partCount_;

    // nets for this hypergraph
    std::vector<std::vector<TIdx>> nets_;

    // nets in which a vector sits
    std::vector<std::vector<TIdx>> netsForVertex_;

    // netDistribution[net] is an array of size partCount that
    // denots how many vertices in a net are in a given part
    std::vector<std::vector<TIdx>> netDistribution_;

    // size of each part
    std::vector<TIdx> partSize_;

    // part of each vertex
    std::vector<TIdx> part_;

    // weight of each vertex
    std::vector<TIdx> weights_;
};

template <typename TIdx = Zee::default_index_type,
          class TMatrix = Zee::DSparseMatrix<>>
class FineGrainHG : public DHypergraph<TIdx> {
  public:
    FineGrainHG(TMatrix& A)
        : DHypergraph<TIdx>(A.nonZeros(), A.getProcs()), A_(A) {
        this->netCount_ = A.getCols() + A.getRows();

        this->nets_.resize(this->netCount_);
        this->netDistribution_.resize(this->netCount_);
        for (auto& distribution : this->netDistribution_)
            distribution.resize(this->partCount_);

        this->netsForVertex_.resize(this->vertexCount_);

        row.resize(this->vertexCount_);
        col.resize(this->vertexCount_);

        TIdx s = 0;
        // we want a fixed ordering of the nonzeros..
        TIdx i = 0;
        for (auto& image : A.getImages()) {
            for (auto& trip : *image) {
                row[i] = trip.row();
                col[i] = trip.col();

                this->part_[i] = s;
                this->partSize_[s]++;

                this->nets_[trip.col()].push_back(i);
                this->netDistribution_[trip.col()][s]++;

                this->nets_[A.getCols() + trip.row()].push_back(i);
                this->netDistribution_[A.getCols() + trip.row()][s]++;

                this->netsForVertex_[i].push_back(trip.col());
                this->netsForVertex_[i].push_back(A.getCols() + trip.row());

                this->weights_[i] = 1;

                i++;
            }

            ++s;
        }
    }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            ZeeLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->partSize_[this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
            this->partSize_[part]++;
        }

        A_.moveNonZero(row[vertex], col[vertex], this->part_[vertex], part);
        this->part_[vertex] = part;
    }

  private:
    TMatrix& A_;
    std::vector<TIdx> row;
    std::vector<TIdx> col;
};

template <typename TIdx = Zee::default_index_type,
          class TMatrix = Zee::DSparseMatrix<>>
class RowNetHG : public DHypergraph<TIdx> {
  public:
    RowNetHG(TMatrix& A) : DHypergraph<TIdx>(A.getCols(), A.getProcs()), A_(A) {
        this->netCount_ = A.getRows();

        this->nets_.resize(this->netCount_);
        this->netDistribution_.resize(this->netCount_);
        for (auto& distribution : this->netDistribution_)
            distribution.resize(this->partCount_);

        this->netsForVertex_.resize(this->vertexCount_);

        TIdx s = 0;
        for (auto& image : A.getImages()) {
            for (auto& trip : *image) {
                this->part_[trip.col()] = s;
                this->weights_[trip.col()]++;
                this->partSize_[s]++;

                this->nets_[trip.row()].push_back(trip.col());
                this->netDistribution_[trip.row()][s]++;

                this->netsForVertex_[trip.col()].push_back(trip.row());
            }

            ++s;
        }
    }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            ZeeLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->partSize_[this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
            this->partSize_[part]++;
        }

        for (TIdx row : this->netsForVertex_[vertex]) {
            A_.moveNonZero(row, vertex, this->part_[vertex], part);
        }
        this->part_[vertex] = part;
    }


  private:
    TMatrix& A_;
};

template <typename TIdx = Zee::default_index_type,
          class TMatrix = Zee::DSparseMatrix<>>
class ColumnNetHG : public DHypergraph<TIdx> {
  public:
    ColumnNetHG(TMatrix& A)
        : DHypergraph<TIdx>(A.getRows(), A.getProcs()), A_(A) {
        this->netCount_ = A.getCols();

        this->nets_.resize(this->netCount_);
        this->netDistribution_.resize(this->netCount_);
        for (auto& distribution : this->netDistribution_)
            distribution.resize(this->partCount_);

        this->netsForVertex_.resize(this->vertexCount_);

        TIdx s = 0;
        // we want a fixed ordering of the nonzeros..
        for (auto& image : A.getImages()) {
            for (auto& trip : *image) {

                this->part_[trip.row()] = s;
                this->weights_[trip.row()]++;
                this->partSize_[s]++;

                this->nets_[trip.col()].push_back(trip.row());
                this->netDistribution_[trip.col()][s]++;

                this->netsForVertex_[trip.row()].push_back(trip.col());
            }

            ++s;
        }
    }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            ZeeLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->partSize_[this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
            this->partSize_[part]++;
        }

        for (TIdx col : this->netsForVertex_[vertex]) {
            A_.moveNonZero(vertex, col, this->part_[vertex], part);
        }
        this->part_[vertex] = part;
    }


  private:
    TMatrix& A_;
};
