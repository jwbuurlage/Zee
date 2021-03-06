#pragma once

#include <algorithm>
#include <random>
#include <vector>

#include "zee.hpp"

namespace Zee {

template <typename TIdx = Zee::default_index_type>
class DHypergraph {
   public:
    DHypergraph(TIdx vertexCount, TIdx partCount)
        : vertexCount_(vertexCount),
          partCount_(partCount),
          partSize_(partCount),
          part_(vertexCount),
          weights_(vertexCount) {}

    virtual ~DHypergraph() = default;

    virtual std::vector<double> partQuality(TIdx v,
                                            std::function<double(TIdx, TIdx)> w,
                                            unsigned int maximumNetSize = 0,
                                            unsigned int netsToConsider = 0) {
        std::vector<double> result(this->partCount_);

        if (maximumNetSize == 0) {
            maximumNetSize = getMaximumNetSize();
        }
        if (netsToConsider == 0) {
            netsToConsider = nets_.size();
        }

        for (auto n : this->netsForVertex_[v]) {
            if (this->nets_[n].size() > maximumNetSize) continue;

            if (this->netsSizeRank_[n] >= netsToConsider) continue;

            for (TIdx i = 0; i < this->partCount_; ++i) {
                result[i] +=
                    w(this->netDistribution_[n][i], this->nets_[n].size());
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
                if (cnt > 0) lambda++;
            result += lambda - 1;
        }

        return result;
    }

    std::vector<std::vector<TIdx>>& getNets() const { return nets_; }

    virtual void reassign(TIdx vertex, TIdx part) = 0;
    virtual void clean() = 0;

    void sortNetsBySize() {
        this->netsSizeRank_.resize(nets_.size());

        std::vector<TIdx> netsBySize(nets_.size());
        std::iota(netsBySize.begin(), netsBySize.end(), 0);

        std::sort(netsBySize.begin(), netsBySize.end(),
                  [this](TIdx lhs, TIdx rhs) {
                      return this->nets_[lhs].size() < this->nets_[rhs].size();
                  });

        for (TIdx i = 0; i < nets_.size(); ++i) {
            this->netsSizeRank_[netsBySize[i]] = i;
        }
    }

    TIdx getMaximumNetSize() {
        if (maximumNetSize_ == 0) {
            for (auto& net : nets_) {
                if (net.size() > maximumNetSize_) maximumNetSize_ = net.size();
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
        storageIndex_.resize(this->vertexCount_);

        TIdx s = 0;
        // we want a fixed ordering of the nonzeros..
        TIdx i = 0;
        for (auto& image : A.getImages()) {
            TIdx k = 0;
            for (auto& trip : *image) {
                row[i] = trip.row();
                col[i] = trip.col();

                this->part_[i] = s;
                this->partSize_[s]++;

                this->storageIndex_[i] = k;

                this->nets_[trip.col()].push_back(i);
                this->netDistribution_[trip.col()][s]++;

                this->nets_[A.getCols() + trip.row()].push_back(i);
                this->netDistribution_[A.getCols() + trip.row()][s]++;

                this->netsForVertex_[i].push_back(trip.col());
                this->netsForVertex_[i].push_back(A.getCols() + trip.row());

                this->weights_[i] = 1;

                ++i;
                ++k;
            }

            ++s;
        }
    }

    void clean() { A_.clean(); }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            JWLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->partSize_[this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
            this->partSize_[part]++;
        }

        storageIndex_[vertex] =
            A_.moveNonZero(storageIndex_[vertex], this->part_[vertex], part);

        this->part_[vertex] = part;
    }

   private:
    TMatrix& A_;
    std::vector<TIdx> row;
    std::vector<TIdx> col;
    // storage index for vertex
    std::vector<TIdx> storageIndex_;
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
        this->storageIndices_.resize(this->vertexCount_);

        TIdx s = 0;
        for (auto& image : A.getImages()) {
            TIdx k = 0;
            for (auto& trip : *image) {
                this->part_[trip.col()] = s;
                this->weights_[trip.col()]++;
                this->partSize_[s]++;

                this->nets_[trip.row()].push_back(trip.col());
                this->netDistribution_[trip.row()][s]++;

                this->netsForVertex_[trip.col()].push_back(trip.row());
                this->storageIndices_[trip.col()].push_back(k);

                ++k;
            }
            ++s;
        }
    }

    void clean() { A_.clean(); }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            JWLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->partSize_[this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
            this->partSize_[part]++;
        }

        for (auto& idx : this->storageIndices_[vertex]) {
            idx = A_.moveNonZero(idx, this->part_[vertex], part);
        }
        this->part_[vertex] = part;
    }

   private:
    TMatrix& A_;
    // storage index for vertex
    std::vector<std::vector<TIdx>> storageIndices_;
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
        this->storageIndices_.resize(this->vertexCount_);

        TIdx s = 0;
        // we want a fixed ordering of the nonzeros..
        for (auto& image : A.getImages()) {
            TIdx k = 0;
            for (auto& trip : *image) {
                this->part_[trip.row()] = s;
                this->weights_[trip.row()]++;
                this->partSize_[s]++;

                this->nets_[trip.col()].push_back(trip.row());
                this->netDistribution_[trip.col()][s]++;

                this->netsForVertex_[trip.row()].push_back(trip.col());
                this->storageIndices_[trip.row()].push_back(k);

                ++k;
            }
            ++s;
        }
    }

    void clean() { A_.clean(); }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            JWLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->partSize_[this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
            this->partSize_[part]++;
        }

        for (auto& idx : this->storageIndices_[vertex]) {
            idx = A_.moveNonZero(idx, this->part_[vertex], part);
        }
        this->part_[vertex] = part;
    }

   private:
    TMatrix& A_;
    // label for every nonzero -- large!
    std::vector<std::vector<TIdx>> storageIndices_;
};

/* A = A_r + A_c
 *
 * B = [ I_n | A_r^T ]
 *     [-------------]
 *     [ A_c | I_m   ]
 *
 * we apply row-net model to this matrix
 * we have n + m vertices
 * we have n + m nets
 * every vertex consists of a number of non-zeros that we
 * should add to it
 */

template <typename TIdx = Zee::default_index_type,
          class TMatrix = Zee::DSparseMatrix<>>
class MediumGrainHG : public DHypergraph<TIdx> {
   public:
    MediumGrainHG(TMatrix& matrix, TIdx procs)
        : DHypergraph<TIdx>(matrix.getCols() + matrix.getRows(), procs),
          matrix_(matrix) {
        // we do the initial partitioning ourselves
        JWAssert(matrix.getProcs() == 1);

        // FIXME move this up
        this->netCount_ = matrix.getCols() + matrix.getRows();
        this->nets_.resize(this->netCount_);
        this->netDistribution_.resize(this->netCount_);
        for (auto& distribution : this->netDistribution_)
            distribution.resize(this->partCount_);
        this->netsForVertex_.resize(this->vertexCount_);

        this->storageIndices_.resize(this->vertexCount_);

        this->split_();

        this->initial_distribution_();
    }

    void clean() { matrix_.clean(); }

    void reassign(TIdx vertex, TIdx part) override {
        if (this->part_[vertex] == part) {
            JWLogError << "Reassigning to own part" << endLog;
            return;
        }

        for (auto net : this->netsForVertex_[vertex]) {
            this->netDistribution_[net][this->part_[vertex]]--;
            this->netDistribution_[net][part]++;
        }

        this->partSize_[this->part_[vertex]] -= this->weights_[vertex];
        this->partSize_[part] += this->weights_[vertex];

        for (auto& idx : this->storageIndices_[vertex]) {
            idx = matrix_.moveNonZero(idx, this->part_[vertex], part);
        }

        this->part_[vertex] = part;
    }

   private:
    void split_() {
        std::vector<TIdx> row_count(matrix_.getRows());
        std::vector<TIdx> col_count(matrix_.getCols());

        for (auto triplet : *matrix_.getImages()[0]) {
            row_count[triplet.row()]++;
            col_count[triplet.col()]++;
        }

        for (TIdx i = 0; i < this->vertexCount_; ++i) {
            this->nets_[i].push_back(i);
            this->part_[i] = 0;
            this->netDistribution_[i][0]++;
        }

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<TIdx> randbool(0, 1);

        TIdx k = 0;
        for (auto triplet : *matrix_.getImages()[0]) {
            TIdx target_vertex = 0;
            TIdx target_net = 0;
            if (row_count[triplet.row()] < col_count[triplet.col()]) {
                target_vertex = matrix_.getCols() + triplet.row();
                target_net = triplet.col();
            } else {
                target_vertex = triplet.col();
                target_net = matrix_.getCols() + triplet.row();
            }
            this->weights_[target_vertex] += 1;
            this->storageIndices_[target_vertex].push_back(k);
            this->nets_[target_net].push_back(target_vertex);
            this->netDistribution_[target_net][0]++;
            this->netsForVertex_[target_vertex].push_back(target_net);
            ++k;
        }

        for (TIdx i = 0; i < this->vertexCount_; ++i) {
            this->partSize_[0] += this->weights_[i];
        }
    }

    void initial_distribution_() {
        using TImage = typename TMatrix::image_type;
        matrix_.getMutableImages().push_back(std::make_shared<TImage>());
        matrix_.setProcs(2);

        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<TIdx> randpart(0, this->partCount_ - 1);

        for (TIdx k = matrix_.getCols(); k < this->vertexCount_; ++k) {
            TIdx target_part = randpart(mt);
            if (target_part != this->part_[k]) this->reassign(k, target_part);
        }
    }

    TMatrix& matrix_;

    // nets_
    // part_
    std::vector<std::vector<TIdx>> storageIndices_;
};

}  // namespace Zee
