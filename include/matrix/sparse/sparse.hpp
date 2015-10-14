/*
This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include <cassert>
#include <cstdint>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <random>
#include <memory>
#include <thread>

#include "../base/base.hpp"
#include "../storage.hpp"
#include "../../common.hpp"
#include "../../logging.hpp"
#include "../../operations/operations.hpp"

namespace Zee {

// Fwd declaring partitioner
template <class TMatrix>
class Partitioner;

template <typename TVal, typename TIdx = uint32_t,
          class CStorage = StorageTriplets<TVal, TIdx>>
class DSparseMatrixImage;

// Matrix Market
namespace matrix_market {
template <typename Derived, typename TVal, typename TIdx>
void load(std::string file, DMatrixBase<Derived, TVal, TIdx>& target);
}

/** A (matrix) triplet is a tuplet \f$(i, j, a_{ij})\f$, representing an entry
 * in a
  * matrix. It is particularly useful in the representation of sparse matrices.
  */
template <typename TVal, typename TIdx = uint32_t>
class Triplet {
  public:
    Triplet(TIdx i, TIdx j, TVal value) {
        i_ = i;
        j_ = j;
        value_ = value;
    };

    ~Triplet() {}

    /** @return the value for this triplet */
    TVal value() const { return value_; }

    /** @return the row position inside the matrix */
    TIdx row() const { return i_; }

    /** @return the column position inside the matrix */
    TIdx col() const { return j_; }

    void setRow(TIdx row)  {
        i_ = row;
    }

    void setCol(TIdx col)  {
        j_ = col;
    }

    void setValue(TVal val)  {
        value_ = val;
    }

  private:
    TIdx i_;
    TIdx j_;
    TVal value_;
};

//-----------------------------------------------------------------------------

/** FIXME DEPRECATED: Different partitioning schemes */
enum class partitioning_scheme { cyclic, block, random, custom };

template <typename Derived, typename TVal, typename TIdx,
          class Image =
              DSparseMatrixImage<TVal, TIdx, StorageTriplets<TVal, TIdx>>>
class DSparseMatrixBase : public DMatrixBase<Derived, TVal, TIdx> {
  public:
    using index_type = TIdx;
    using value_type = TVal;
    using image_type = Image;

    using Base = DMatrixBase<Derived, TVal, TIdx>;
    using Base::operator=;

    DSparseMatrixBase(TIdx rows, TIdx cols) : Base(rows, cols) {}

    DSparseMatrixBase(std::string file, TIdx procs = 0) : Base(0, 0) {
        setDistributionScheme(partitioning_scheme::cyclic, procs);
        matrix_market::load(file, *(Derived*)this);
    }

    /** Construct a matrix from a set of triplets */
    template <typename TInputIterator>
    void setFromTriplets(const TInputIterator& begin,
                         const TInputIterator& end) {
        subs_.clear();

        nz_ = 0;

        ZeeAssert(this->getProcs() > 0);

        if (partitioning_ == partitioning_scheme::custom) {
            if (!distributionLambda_) {
                ZeeLogError
                    << "Trying to apply a custom partitioning, but"
                       " no distribution function was set. The matrix remains"
                       " uninitialized." << endLog;
                return;
            }
        }

        for (TIdx i = 0; i < this->getProcs(); ++i)
            subs_.push_back(std::make_shared<Image>());

        // FIXME: only if partitioning is random
        std::random_device rd;
        std::mt19937 mt(rd());
        std::uniform_int_distribution<TIdx> randproc(0, this->getProcs() - 1);

        // FIXME change order of switch and for
        for (TInputIterator it = begin; it != end; it++) {
            TIdx target_proc = 0;
            switch (partitioning_) {
            case partitioning_scheme::cyclic:
                target_proc = (*it).row() % this->getProcs();
                break;

            case partitioning_scheme::block:
                target_proc =
                    (this->getProcs() * (*it).row()) / this->getRows();
                break;

            case partitioning_scheme::random:
                target_proc = randproc(mt);
                break;

            case partitioning_scheme::custom:
                target_proc = distributionLambda_((*it).row(), (*it).col());
                break;

            default:
                // Fall back to 1D cyclic
                target_proc = (*it).row() % this->getProcs();
                break;
            }

            subs_[target_proc]->pushTriplet(*it);
            nz_++;
        }

        initialized_ = true;
    }

    /** Sets the distribution scheme for this matrix */
    void setDistributionScheme(partitioning_scheme partitioning, TIdx procs) {
        this->partitioning_ = partitioning;
        this->procs_ = procs;
    }

    /** Sets the distribution scheme for this matrix. The function should be of
     * the form
      * \f[ f: Z_m \times Z_n \to Z_p \f]
      * where m is the number of columns, n is the number of rows, and p is the
     * number
      * of processors of this matrix */
    void setDistributionFunction(
        std::function<TIdx(TIdx, TIdx)> distributionLambda) {
        distributionLambda_ = distributionLambda;
    }

    /** @return the number of non-zero entries in the matrix */
    TIdx nonZeros() const { return nz_; }

    /** Obtain a list of images */
    const std::vector<std::shared_ptr<Image>>& getImages() const {
        return subs_;
    }

    std::vector<std::shared_ptr<Image>>& getMutableImages() { return subs_; }

    // this is kind of like a reduce in mapreduce, implementing this such that
    // we can get some sample code going
    // perhaps think about pregel-like approach as well
    // An alternative function signature could be:
    //   template<typename TReturn, typename TFunc>
    //   std::vector<TReturn> compute(TFunc func) const
    // which would construct a separate function for each lambda, but remove
    // overhead
    // of using std::function.
    template <typename TReturn>
    std::vector<TReturn>
    compute(std::function<TReturn(std::shared_ptr<image_type>)> func) const {
        auto result = std::vector<TReturn>(this->getProcs());

        // capture function and result by-reference
        auto push_result = [&func, &result](
            TIdx proc, const std::shared_ptr<image_type>& image) mutable {
            result[proc] = func(image);
        };

        std::vector<std::thread> threads(this->subs_.size());

        TIdx p = 0;
        for (const auto& image : this->subs_) {
            threads[p] = std::thread(push_result, p, image);
            ++p;
        }

        for (auto& t : threads) {
            t.join();
        }

        return result;
    }

    // template specialization for void which does not return anything
    void compute(
        std::function<void(std::shared_ptr<image_type>, TIdx s)> func) const {
        // capture function and result by-reference
        std::vector<std::thread> threads;

        TIdx s = 0;
        for (auto& image : this->subs_) {
            threads.push_back(std::thread(func, image, s++));
        }

        for (auto& t : threads) {
            t.join();
        }
    }

    /** Returns the load imbalance of the current partitioning.
     * Load imbalance is defined as:
     * \f[ \tilde{\epsilon} = \max_{i \in P} \frac{p \cdot |A_i|}{|A|} \f]
     * and should be smaller than some predetermined value \f$\epsilon\f$.
     */
    virtual double loadImbalance() const {
        double eps = 1.0;
        for (auto& pimg : this->subs_) {
            double eps_i = ((double)this->getProcs() * pimg->nonZeros()) /
                           this->nonZeros();
            if (eps_i > eps) {
                eps = eps_i;
            }
        }

        return eps;
    }

    /** Returns the communication volume of the current Partitioning
      * Let \f$\lambda_i\f$ be the number of processors in a (non-empty)
      * row \f$i\f$, and \f$\mu_j\f$ be the number of processors in a
      * (non-empty) column \f$j\f$. Then the total communication volume is:
      * \f[ V = \sum_i (\lambda_i - 1) + \sum_j (\mu_j - 1) \f]
      * */
    virtual TIdx communicationVolume() const {
        // here we assume that v_i is owned by a processor holding
        // a_ik =/= 0 for some k, and u_j is owned by a processor
        // holding a_kj =/= 0 for some k.
        //
        // We then ask the images themselves how much communciation
        // is necessary for spmv communication.
        //
        // We preprocess a 'rows' and 'cols' vector in each image.
        // e.g. for i in rows, send (1, i) to (i % p) proc.
        //
        // proc adds 1 to lambda_i for each O(n / p) i's it owns.
        //
        // In the end lambda's are distributed over procs, let them compute
        // partial sums and send it upwards.

        std::atomic<TIdx> V{0};

        // {lambda,mu}_i = {lambda,mu}[i % p][i / p];
        std::vector<std::vector<atomic_wrapper<TIdx>>> lambda(
            this->getProcs(), std::vector<atomic_wrapper<TIdx>>(
                                  this->getRows() / this->getProcs() + 1));
        std::vector<std::vector<atomic_wrapper<TIdx>>> mu(
            this->getProcs(), std::vector<atomic_wrapper<TIdx>>(
                                  this->getCols() / this->getProcs() + 1));

        // FIXME: parallelize, using generalized compute
        TIdx s = 0;
        for (auto& pimg : this->subs_) {
            for (auto key_count : pimg->getRowSet()) {
                auto i = key_count.first;
                lambda[i % this->getProcs()][i / this->getProcs()].a += 1;
            }

            for (auto key_count : pimg->getColSet()) {
                auto j = key_count.first;
                mu[j % this->getProcs()][j / this->getProcs()].a += 1;
            }

            ++s;
        }

        // now let each proc compute partial sum

        // sum over lambdas
        // FIXME: parallelize, using generalized compute
        // FIXME: can also parallelize.. O(n) -> O(n / p)
        // with extra superstep
        std::vector<TIdx> lambda_s(this->getProcs(), 0);
        for (TIdx proc = 0; proc < this->getProcs(); ++proc) {
            TIdx V_s = 0;

            for (TIdx i = 0; i < lambda[proc].size(); ++i) {
                if (lambda[proc][i].a > 1) {
                    V_s += lambda[proc][i].a - 1;
                }
            }

            for (TIdx i = 0; i < mu[proc].size(); ++i) {
                if (mu[proc][i].a > 1) {
                    V_s += mu[proc][i].a - 1;
                }
            }

            V += V_s;
        }

        return V;
    }

    void resetImages(std::vector<std::unique_ptr<Image>>& new_images) {
        // update the number of processors
        this->procs_ = new_images.size();
        this->subs_.resize(this->getProcs());

        for (TIdx i = 0; i < new_images.size(); ++i)
            this->subs_[i].reset(new_images[i].release());

        // update nz_
        this->nz_ = 0;
        for (auto& pimg : this->subs_) {
            this->nz_ += pimg->nonZeros();
        }

        this->initialized_ = true;
    }

    /** Obtain the number of nonzeros in a column */
    TIdx getColumnWeight(TIdx j) const {
        // TODO precompute
        // TODO optimize
        auto columnCount = [j](std::shared_ptr<image_type> image) -> TIdx {
            TIdx count = 0;
            for (auto& trip : *image)
                if (trip.col() == j)
                    count++;
            return count;
        };

        auto counts = this->template compute<TIdx>(columnCount);
        TIdx count = std::accumulate(counts.begin(), counts.end(), (TIdx)0);

        return count;
    }

    /** Obtain a spy image of the sparse matrix */
    void spy(std::string title = "anonymous", bool show = false) {
        using std::endl;

        std::stringstream ss;
        ss << "data/spies/" << title << ".mtx";
        auto filename = ss.str();
        int i = 1;
        while (fileExists(filename)) {
            ss.str("");
            ss.clear();
            ss << "data/spies/" << title << "_" << i++ << ".mtx";
            filename = ss.str();
        }
        std::ofstream fout(filename);

        fout << "%%MatrixMarket matrix coordinate integer general" << endl;

        fout << "% Matrix sparsity:      " << std::fixed << std::setprecision(4)
             << this->nonZeros() / static_cast<double>(this->size()) << endl;

        fout << "% Load imbalance:       " << std::fixed << std::setprecision(4)
             << this->loadImbalance() << endl;

        fout << "% Communication Volume: " << this->communicationVolume()
             << endl;
        fout << title << endl;

        fout << this->getRows() << " " << this->getCols() << " "
             << this->nonZeros() << endl;

        TIdx s = 0;
        for (auto& image : this->subs_) {
            for (auto& triplet : *image) {
                fout << triplet.row() << " " << triplet.col() << " " << s
                     << endl;
            }
            ++s;
        }

        ZeeLogInfo << "Spy saved to file: " << filename << Logger::end();

        if (show) {
            auto command = "./script/plot.py --showfile " + filename;
            std::system(command.c_str());
        }
    }

    bool isInitialized() const { return this->initialized_; }

  protected:
    TIdx nz_;
    partitioning_scheme partitioning_;
    std::vector<std::shared_ptr<Image>> subs_;
    std::function<TIdx(TIdx, TIdx)> distributionLambda_;
    bool initialized_ = false;
};

/** The class DSparseMatrix is a distributed matrix type inspired by
  * Eigen's SparseMatrix. It is distributed over multiple processing units,
  */
template <
    typename TVal = default_scalar_type, typename TIdx = default_index_type,
    class Image = DSparseMatrixImage<TVal, TIdx, StorageTriplets<TVal, TIdx>>>
class DSparseMatrix : public DSparseMatrixBase<DSparseMatrix<TVal, TIdx, Image>,
                                               TVal, TIdx, Image> {
    using Base =
        DSparseMatrixBase<DSparseMatrix<TVal, TIdx, Image>, TVal, TIdx, Image>;
    using Base::operator=;

  public:
    using index_type = typename Base::index_type;
    using value_type = typename Base::value_type;
    using image_type = typename Base::image_type;

    /** Initialize from .mtx format */
    DSparseMatrix(std::string file, TIdx procs = 1) : Base(file, procs) {}

    /** Initialize an (empty) sparse (rows x cols) matrix */
    DSparseMatrix(TIdx rows, TIdx cols, TIdx procs = 0) : Base(rows, cols) {
        this->setDistributionScheme(partitioning_scheme::cyclic, procs);
    }

    DSparseMatrix() : DSparseMatrix<TVal, TIdx>(0, 0) {}

    /** Default deconstructor */
    ~DSparseMatrix() {}

    /** Move constructor */
    DSparseMatrix(DSparseMatrix&& o) = default;

};

// Owned by a processor. It is a submatrix, which holds the actual
// data, the 'global' DSparseMatrix can be seen as the sum of these images.
//
// BRAINSTORM:
// - What if every matrix image had a 'real', physical chunk of memory and
//   a compute unit, whatever that would be
// - Then we could let programs execute code on these images (as in, lambdas)
// - Furthermore, different compute units could execute code remotely without
//   the 'remote code' explicitely knowing about this.
// - Is this flexible enough, and would this work on e.g. Cartesius?
// - I want to try and implement a version of this

template <typename TVal, typename TIdx, class CStorage>
class DSparseMatrixImage {
  public:
    /** Default constructor */
    DSparseMatrixImage() : storage_(new CStorage()) {
        // FIXME: Switch cases between different storage types
    }

    // Because we are using a unique pointer we need to move ownership
    // upon copying
    DSparseMatrixImage(DSparseMatrixImage&& other)
        : storage_(std::move(other.storage_)) {}

    ~DSparseMatrixImage() = default;

    // FIXME rename to erase
    void popElement(TIdx element) {
        auto t = storage_->popElement(element);
        rowset_.lower(t.row());
        colset_.lower(t.col());
    }

    // rename to push
    void pushTriplet(Triplet<TVal, TIdx> t) {
        if (!storage_) {
            ZeeLogError << "Can not push triplet without storage." << endLog;
            return;
        }
        rowset_.raise(t.row());
        colset_.raise(t.col());
        storage_->pushTriplet(t);
    }

    void setLocalIndices(std::vector<TIdx>&& localIndicesV,
                         std::vector<TIdx>&& localIndicesU) {
        localIndicesV_ = localIndicesV;
        localIndicesU_ = localIndicesU;

        // obtain number of local components
        numLocalV_ = localIndicesV.size();
        numLocalU_ = localIndicesU.size();

        // add necessary remote indices to list of local indices
        computeLocalIndices_(localIndicesV_, colset_);
        computeLocalIndices_(localIndicesU_, rowset_);
    }

    void localizeStorage() {
        ZeeLogDebug << "Localize storage" << endLog;
        // 1. construct inverse map
        std::map<TIdx, TIdx> globalToLocalV;
        std::map<TIdx, TIdx> globalToLocalU;
        for (TIdx idx = 0; idx < localIndicesV_.size(); ++idx) {
            globalToLocalV[localIndicesV_[idx]] = idx;
        }
        for (TIdx idx = 0; idx < localIndicesU_.size(); ++idx) {
            globalToLocalU[localIndicesU_[idx]] = idx;
        }
        ZeeLogVar(globalToLocalV);
        ZeeLogVar(globalToLocalU);

        // 2. let storage localize itself
        storage_->localize(globalToLocalV, globalToLocalU);

        localizedStorage_ = true;
    }

    const counted_set<TIdx>& getRowSet() const { return rowset_; }

    const counted_set<TIdx>& getColSet() const { return colset_; }

    using iterator = typename CStorage::it_traits::iterator;

    iterator begin() const { return storage_->begin(); }

    iterator end() const { return storage_->end(); }

    /** @return The number of nonzeros in this image */
    TIdx nonZeros() const { return storage_->size(); }

    /** @return The \f$i\f$-th element in this image */
    Triplet<TVal, TIdx> getElement(TIdx i) const {
        return storage_->getElement(i);
    }

    std::vector<TIdx>& getLocalIndicesU() { return localIndicesU_; }
    std::vector<TIdx>& getLocalIndicesV() { return localIndicesV_; }

    TIdx getNumLocalU() { return numLocalU_; }
    TIdx getNumLocalV() { return numLocalV_; }

    std::vector<TIdx>& getRemoteOwnersU() { return remoteOwnersU_; }
    std::vector<TIdx>& getRemoteOwnersV() { return remoteOwnersV_; }

    bool localizedStorage() const { return localizedStorage_; }

  private:
    void computeLocalIndices_(std::vector<TIdx>& partialLocalIndices,
            const counted_set<TIdx>& countedSetOfIndices) {
        TIdx numLocal = partialLocalIndices.size();
        TIdx localIdx = 0;
        for (auto& idx : countedSetOfIndices) {
            while (localIdx < numLocal &&
                   partialLocalIndices[localIdx] < idx.first) {
                // we do not own anything in the column
                // localIndicesV[localIdxV], but hold the index because of the
                // other vector, we *may still have idx.first*, so we increase
                // the local index until we reach idx.first
                localIdx++;
            }

            if (localIdx == numLocal ||
                idx.first != partialLocalIndices[localIdx]) {
                // we dont own it
                partialLocalIndices.push_back(idx.first);
            }
            else {
                // we own it
                localIdx++;
            }
        }
    }

    /** We delegate the storage to a superclass (to simplify choosing
     * a storage mechanism) */
    std::unique_ptr<CStorage> storage_;

    // local vector indices
    std::vector<TIdx> localIndicesU_;
    std::vector<TIdx> localIndicesV_;
    // number of local components U
    TIdx numLocalU_;
    TIdx numLocalV_;
    // where to obtain missing non-local components
    std::vector<TIdx> remoteOwnersU_;
    std::vector<TIdx> remoteOwnersV_;

    /** A set (with counts) that stores the non-empty rows in this image */
    counted_set<TIdx> rowset_;
    /** A set (with counts) that stores the non-empty columns in this image */
    counted_set<TIdx> colset_;

    // whether we already localized storage
    bool localizedStorage_ = false;

    /** We hold a reference to the other images */
    // FIXME implement and use
    std::vector<std::weak_ptr<DSparseMatrixImage<TVal, TIdx, CStorage>>>
        images_;
};


//-----------------------------------------------------------------------------
// LOGGING

template <typename TVal, typename TIdx>
Logger& operator<<(Logger& lhs, const Triplet<TVal, TIdx>& rhs) {
    lhs << "{" << rhs.row() << ", " << rhs.col() << ", " << rhs.value() << "}";
    return lhs;
}

} // namespace Zee
