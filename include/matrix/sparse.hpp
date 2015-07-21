/*
File: include/sparse_matrix.h

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

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <algorithm>
#include <vector>
#include <random>
#include <memory>
#include <thread>

#include "matrix/base.hpp"
#include "matrix/dense.hpp"
#include "matrix/storage.hpp"
#include "color_output.hpp"
#include "logging.hpp"
#include "common.hpp"

namespace Zee {

// Fwd declaring partitioner
template <class TMatrix>
class Partitioner;

template <typename TVal, typename TIdx = uint32_t,
         class CStorage = StorageTriplets<TVal, TIdx>>
class DSparseMatrixImage;

/** A (matrix) triplet is a tuplet \f$(i, j, a_{ij})\f$, representing an entry in a
  * matrix. It is particularly useful in the representation of sparse matrices.
  */
template <typename TVal, typename TIdx = uint32_t>
class Triplet
{
    public:
        Triplet(TIdx i, TIdx j, TVal value)
        {
            _i = i;
            _j = j;
            _value = value;
        };

        ~Triplet() { }

        /** @return the value for this triplet */
        inline TVal value() const { return _value; }

        /** @return the row position inside the matrix */
        inline TIdx row() const { return _i; }

        /** @return the column position inside the matrix */
        inline TIdx col() const { return _j; }

    private:
        TIdx _i;
        TIdx _j;
        TVal _value;
};

//-----------------------------------------------------------------------------

/** Different partitioning schemes */
enum class Partitioning
{
    cyclic,
    block,
    random,
    custom
};

/** The class DSparseMatrix is a distributed matrix type inspired by
  * Eigen's SparseMatrix. It is distributed over multiple processing units,
  */
template <typename TVal, typename TIdx = uint32_t,
         class Image = DSparseMatrixImage<TVal, TIdx,
            StorageTriplets<TVal, TIdx>>>
class DSparseMatrix : public DMatrixBase<TVal, TIdx>
{
    public:
        using image_type = Image;
        using index_type = TIdx;
        using value_type = TVal;

        /** Initialize an (empty) sparse (rows x cols) oatrix */
        DSparseMatrix(TIdx rows, TIdx cols) :
            DMatrixBase<TVal, TIdx>(rows, cols)
        {
            _partitioning = Partitioning::cyclic;
        }

        /** Default deconstructor */
        ~DSparseMatrix() { }

        /** Move constructor */
        DSparseMatrix(DSparseMatrix&& o) = default;

        /** Sets the distribution scheme for this matrix */
        void setDistributionScheme(Partitioning partitioning,
                TIdx procs)
        {
            _partitioning = partitioning;
            this->_procs = procs;
        }

        /** Multiply a sparse matrix with a dense vector */
        DVector<TVal, TIdx> operator*(const DVector<TVal, TIdx>& v) const
        {
            DVector<TVal, TIdx> u(this->rows());
            // how is u distributed, depends on left hand side.. very weird
            // construction, but interesting. Should return expression template
            // here
            // distributed shit
            return u;
        }

        /** @return the number of non-zero entries in the matrix */
        TIdx nonZeros() const
        {
            return _nz;
        }

        // this is kind of like a reduce in mapreduce, implementing this such that
        // we can get some sample code going
        // perhaps think about pregel-like approach as well
        // An alternative function signature could be:
        //   template<typename TReturn, typename TFunc>
        //   std::vector<TReturn> compute(TFunc func) const
        // which would construct a separate function for each lambda, but remove overhead
        // of using std::function.
        template<typename TReturn>
        std::vector<TReturn> compute(std::function<TReturn(std::shared_ptr<image_type>)> func) const
        {
            auto result = std::vector<TReturn>(this->_procs);

            // capture function and result by-reference
            auto push_result = [&func, &result] (TIdx proc,
                    const std::shared_ptr<image_type>& image) mutable {
                    result[proc] = func(image);
                };

            std::vector<std::thread> threads(this->_procs);

            TIdx p = 0;
            for (const auto& image : _subs) {
                threads[p] = std::thread(push_result, p, image);
                ++p;
            }

            for (auto& t : threads) {
                t.join();
            }

            return result;
        }

        /** Returns the load imbalance of the current partitioning.
         * Load imbalance is defined as:
         * \f[ \tilde{\epsilon} = \max_{i \in P} \frac{p \cdot |A_i|}{|A|} \f]
         * and should be smaller than some predetermined value \f$\epsilon\f$.
         */
        double loadImbalance()
        {
            double eps = 1.0;
            for (auto& pimg : _subs)
            {
                double eps_i = ((double)this->_procs * pimg->nonZeros())
                    / this->nonZeros();
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
        TIdx communicationVolume()
        {
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

            atomic<TIdx> V { 0 };

            // {lambda,mu}_i = {lambda,mu}[i % p][i / p];
            std::vector<std::vector<atomic_wrapper<TIdx>>> lambda(this->_procs,
                    std::vector<atomic_wrapper<TIdx>>(this->_rows / this->_procs + 1));
            std::vector<std::vector<atomic_wrapper<TIdx>>> mu(this->_procs,
                    std::vector<atomic_wrapper<TIdx>>(this->_cols / this->_procs + 1));

            // FIXME: parallelize, using generalized compute
            TIdx s = 0;
            for (auto& pimg : _subs) {
                for (auto key_count : pimg->getRowSet()) {
                    auto i = key_count.first;
                    lambda[i % this->_procs][i / this->_procs]._a += 1;
                }

                for (auto key_count : pimg->getColSet()) {
                    auto j = key_count.first;
                    mu[j % this->_procs][j / this->_procs]._a += 1;
                }

                ++s;
            }

            // now let each proc compute partial sum

            // sum over lambdas
            // FIXME: parallelize, using generalized compute
            // FIXME: can also parallelize.. O(n) -> O(n / p)
            // with extra superstep
            std::vector<int> lambda_s(this->_procs, 0);
            for (TIdx proc = 0; proc < this->_procs; ++proc) {
                TIdx V_s = 0;

                for (int i = 0; i < lambda[proc].size(); ++i) {
                    if (lambda[proc][i]._a > 1) {
                        V_s += lambda[proc][i]._a - 1;
                    }
                }

                for (int i = 0; i < mu[proc].size(); ++i) {
                    if (mu[proc][i]._a > 1) {
                        V_s += mu[proc][i]._a - 1;
                    }
                }

                V += V_s;
            }

            return V;
        }

        /** Construct a matrix from a set of triplets */
        template<typename TInputIterator>
        void setFromTriplets(
            const TInputIterator& begin,
            const TInputIterator& end)
        {
            _subs.clear();

            _nz = 0;

            for (TIdx i = 0; i < this->procs(); ++i)
                _subs.push_back(std::make_shared<Image>());

            // FIXME: only if partitioning is random
            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<TIdx> randproc(0, this->procs() - 1);

            for (TInputIterator it = begin; it != end; it++)
            {
                TIdx target_proc = 0;
                switch (_partitioning)
                {
                case Partitioning::cyclic:
                    target_proc = (*it).row() % this->procs();
                    break;

                case Partitioning::block:
                    target_proc = (this->procs() * (*it).row()) / this->rows();
                    break;

                case Partitioning::random:
                    target_proc = randproc(mt);
                    break;

                default:
                    // Fall back to 1D cyclic
                    target_proc = (*it).row() % this->procs();
                    break;
                }

                _subs[target_proc]->pushTriplet(*it);
                _nz++;
            }
        }

        void resetImages(std::vector<std::unique_ptr<Image>>& new_images)
        {
            // update the number of processors
            this->_procs = new_images.size();
            _subs.resize(this->_procs);

            for(TIdx i = 0; i < new_images.size(); ++i)
                _subs[i].reset(new_images[i].release());

            // update _nz
            _nz = 0;
            for (auto& pimg : _subs) {
                _nz += pimg->nonZeros();
            }
        }

        /** Obtain the number of nonzeros in a column */
        int getColumnWeight(TIdx j) const
        {
            // TODO precompute
            // TODO optimize
            auto columnCount =
                [j] (std::shared_ptr<image_type> image) -> TIdx {
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

        /** Obtain a list of images */
        const std::vector<std::shared_ptr<Image>>& getImages() const
        {
            return _subs;
        }

        std::vector<std::shared_ptr<Image>>& getMutableImages()
        {
            return _subs;
        }

        void spy(std::string title = "anonymous")
        {
            using std::endl;

            std::stringstream ss;
            ss << "data/spies/" << title << ".emtx";
            auto filename = ss.str();
            int i = 1;
            while(fileExists(filename)) {
                ss.str("");
                ss.clear();
                ss << "data/spies/" << title << "_" << i++ << ".emtx";
                filename = ss.str();
            }
            std::ofstream fout(filename);

            fout << "%%Extended-MatrixMarket distributed-matrix coordinate pattern general" << endl;

            fout << "% Matrix sparsity:      " <<
                std::fixed << std::setprecision(4) <<
                this->nonZeros() / static_cast<double>(this->size()) << endl;

            fout << "% Load imbalance:       " <<
                std::fixed << std::setprecision(4) <<
                this->loadImbalance() << endl;

            fout << "% Communication Volume: " <<
                this->communicationVolume() << endl;
            fout << title << endl;

            fout << this->rows() << " " << this->cols() << " " << this->nonZeros() <<
                " " << this->procs() << endl;

            TIdx s = 0;
            for (auto& image : _subs) {
                for (auto& triplet : *image) {
                    fout << triplet.row() << " " << triplet.col() << " " << s << endl;
                }
                ++s;
            }

            ZeeLogInfo << "Spy saved to file: " << filename << Logger::end();
        }

    private:
        TIdx _nz;
        Partitioning _partitioning;
        std::vector<std::shared_ptr<Image>> _subs;
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
class DSparseMatrixImage
{
    public:

        /** Default constructor */
        DSparseMatrixImage() :
            _storage(new CStorage())
        {
            // FIXME: Switch cases between different storage types
        }

        // Because we are using a unique pointer we need to move ownership
        // upon copying
        DSparseMatrixImage(DSparseMatrixImage&& other) :
            _storage(std::move(other._storage)) { }

        ~DSparseMatrixImage() = default;

        void popElement(TIdx element) {
            auto t = _storage->popElement(element);
            _rowset.lower(t.row());
            _colset.lower(t.col());
        }

        void pushTriplet(Triplet<TVal, TIdx> t) {
            assert(_storage);
            _rowset.raise(t.row());
            _colset.raise(t.col());
            _storage->pushTriplet(t);
        }

        const counted_set<TIdx>& getRowSet() const {
            return _rowset;
        }

        const counted_set<TIdx>& getColSet() const {
            return _colset;
        }

        using iterator = typename CStorage::it_traits::iterator;

        iterator begin() const
        {
            return _storage->begin();
        }

        iterator end() const
        {
            return _storage->end();
        }

        /** @return The number of nonzeros in this image */
        TIdx nonZeros() const
        {
            return _storage->size();
        }

        /** @return The \f$i\f$-th element in this image */
        Triplet<TVal, TIdx> getElement(TIdx i) const {
            return _storage->getElement(i);
        }

    private:
        /** We delegate the storage to a superclass (to simplify choosing
         * a storage mechanism) */
        std::unique_ptr<CStorage> _storage;

        /** A set (with counts) that stores the non-empty rows in this image */
        counted_set<TIdx> _rowset;
        /** A set (with counts) that stores the non-empty columns in this image */
        counted_set<TIdx> _colset;

        /** We hold a reference to the other images */
        // FIXME implement and use
        std::vector<std::weak_ptr<DSparseMatrixImage<TVal, TIdx, CStorage>>> _images;
};

//-----------------------------------------------------------------------------
// Convenience functions (MATLAB syntax)
//-----------------------------------------------------------------------------

/** Create identity matrix as sparse matrix */
template <typename TIdx>
DSparseMatrix<double> eye(TIdx n, TIdx procs)
{
    std::vector<Triplet<double, TIdx>> coefficients;
    coefficients.reserve(n);

    for (TIdx i = 0; i < n; ++i)
        coefficients.push_back(Triplet<double, TIdx>(i, i, 1.0));

    DSparseMatrix<double, TIdx> A(n, n);
    A.setDistributionScheme(Partitioning::cyclic, procs);

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

/** Create a random sparse (n x m) matrix */
template <typename TIdx>
DSparseMatrix<double, TIdx> rand(TIdx m, TIdx n, TIdx procs, double density)
{
    std::vector<Triplet<double, TIdx>> coefficients;
    coefficients.reserve((int)(n * m * density));

    double mu = 1.0 / density + 0.5;
    double sigma = 0.5 * mu;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::normal_distribution<double> gauss(mu, sigma);

    TIdx j = static_cast<int>(gauss(mt)) / 2;
    j = std::max(j, (TIdx)1);
    TIdx i = 0;
    while (true)
    {
        coefficients.push_back(Triplet<double, TIdx>(i, j,
                (1.0 + 10.0 * dist(mt))));

        int offset = static_cast<int>(gauss(mt));
        offset = std::max(offset, 1);
        j += offset;

        while (j >= n)
        {
            j = j - n;
            i++;
        }

        if (i >= m)
            break;
    }

    DSparseMatrix<double, TIdx> A(m, n);
    A.setDistributionScheme(Partitioning::random, procs);

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

//-----------------------------------------------------------------------------

} // namespace Zee
