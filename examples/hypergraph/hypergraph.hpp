#include <zee.hpp>

template <typename TIdx = Zee::default_index_type>
class DHypergraph {
  public:
    DHypergraph(TIdx vertexCount, TIdx partCount)
        : vertexCount_(vertexCount), partCount_(partCount),
          partSize_(partCount), part_(vertexCount) {}

    virtual std::vector<TIdx> getNeighbourDistribution(TIdx v) = 0;

    std::vector<std::vector<TIdx>>& getNets() const { return nets_; }
    virtual void reassign(TIdx vertex, TIdx part) = 0;

    TIdx getVertexCount() const { return vertexCount_; }

    TIdx getPart(TIdx v) const { return part_[v]; }
    TIdx getPartSize(TIdx p) const { return partSize_[p]; }

  protected:
    // number of vertices in this hypergraph
    TIdx vertexCount_;

    // number of nets in this hypergraph
    TIdx netCount_;

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
};

template <typename TIdx = Zee::default_index_type, class TMatrix = Zee::DSparseMatrix<>>
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

    std::vector<TIdx> getNeighbourDistribution(TIdx v) override {
        std::vector<TIdx> result(this->partCount_);
        for (auto n : this->netsForVertex_[v]) {
            for (auto u : this->nets_[n]) {
                result[this->part_[u]]++;
            }
        }
        return result;
    }

  private:
    TMatrix& A_;
    std::vector<TIdx> row;
    std::vector<TIdx> col;
};

template <typename TIdx = Zee::default_index_type>
class RowNetHG : DHypergraph<TIdx> {};

template <typename TIdx = Zee::default_index_type>
class ColumnNetHG : DHypergraph<TIdx> {};
