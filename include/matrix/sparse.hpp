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
#include <algorithm>
#include <vector>
#include <random>

#include "matrix/base.hpp"
#include "matrix/dense.hpp"
#include "matrix/storage.hpp"
#include "color_output.hpp"

namespace Zee {

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::vector;

template <typename TVal, typename TIdx = uint32_t,
         class Storage = StorageTriplets<TVal, TIdx>>
class DSparseMatrixImage;

/** A (matrix) triplet is a tuplet (i, j, a_ij), representing an entry in a
  * matrix. It is particularly useful in the representation of sparse matrices.
  */
template <typename TVal, typename TIdx = int32_t>
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
            
        /** Initialize an (empty) sparse (rows x cols) oatrix */
        DSparseMatrix(TIdx rows, TIdx cols) :
            DMatrixBase<TVal, TIdx>(rows, cols)
        {
            _partitioning = Partitioning::cyclic;
        }

        /** Default deconstructor */
        ~DSparseMatrix() = default;

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

        /** Returns the load imbalance of the current partitioning */
        double loadImbalance()
        {
            return -1.0;
        }

        /** Returns the communication volume of the current partitioning */
        TIdx communicationVolume()
        {
            return -1;
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
                _subs.push_back(std::make_unique<Image>());

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

                default:
                    // Fall back to 1D cyclic
                    target_proc = (*it).row() % this->procs();
                    break;
                }

                _subs[target_proc]->pushTriplet(*it);
                _nz++;
            }
        }

        /** Obtain a list of images */
        const vector<unique_ptr<Image>>& getImages() const 
        {
            return _subs;
        }

        void spy()
        {
            cout << endl << "Sparsity pattern of matrix" << endl;
            cout << " - rows: " << this->rows() << endl;
            cout << " - cols: " << this->cols() << endl;
            cout << " - non-zeros: " << this->nonZeros() << endl;
            cout << " - Matrix sparsity: " <<
                this->nonZeros() / static_cast<double>(this->size()) << endl;

            TIdx max_size = 100;
            if (this->rows() > max_size || this->cols() > max_size)
            {
                // FIXME: should just down scale large matrix to
                // 'any point in region'
                cout << "Can not spy large matrices" << endl;
                return;
            }

            if (this->procs() > 4)
            {
                cout << "Spy only supports showing <= 4 procs" << endl;
                return;
            }

            vector<int> output(this->size(), -1);
            for (TIdx i = 0; i < this->size(); ++i)
                output[i] = -1;

            int p = 0;
            for (auto& image : _subs)
            {
                for (auto& triplet : *image)
                    output[triplet.row() * this->cols() + triplet.col()] = p;
                ++p;
            }

            for (TIdx i = 0; i < this->cols() + 2; ++i)
                cout << "__";
            cout << endl;

            static const Color colors[4] = {
                Color::red,
                Color::blue,
                Color::yellow,
                Color::green
            };

            for (TIdx i = 0; i < this->rows(); ++i) {
                cout << "| ";
                for (TIdx j = 0; j < this->cols(); ++j)
                    if (output[i * this->cols() + j] >= 0) {
                        cout << colorOutput(
                                colors[output[i * this->cols() + j]]);
                        cout << output[i * this->cols() + j] << ' ';
                        cout << colorOutput(Color::clear);
                    }
                    else
                        cout << "  ";
                cout << " |" << endl;
            }

            for (TIdx i = 0; i < this->cols() + 2; ++i)
                cout << "--";
            cout << endl;
        }

    private:
        TIdx _nz;

        Partitioning _partitioning;
        vector<unique_ptr<Image>> _subs;
};

// Owned by a processor. It is a submatrix, which holds the actual
// data, the 'global' DSparseMatrix can be seen as the sum of these images.
template <typename TVal, typename TIdx, class Storage>
class DSparseMatrixImage
{
    public:
        /** Default constructor */
        DSparseMatrixImage() :
            _storage(new Storage())
        {
            // FIXME: Switch cases between different storage types
        }

        // Because we are using a unique pointer we need to move ownership
        // upon copying
        DSparseMatrixImage(DSparseMatrixImage&& other) :
            _storage(std::move(other._storage)) { }

        ~DSparseMatrixImage() = default;

        void pushTriplet(Triplet<TVal, TIdx> t) {
            assert(_storage);
            _storage->pushTriplet(t);
        }

        using iterator = typename Storage::it_traits::iterator;

        iterator begin() const
        {
            return _storage->begin();
        }

        iterator end() const
        {
            return _storage->end();
        }

    private:
        // Triplets, CRS or CCS
        unique_ptr<Storage> _storage;
};

//-----------------------------------------------------------------------------
// Convenience functions (MATLAB syntax)
//-----------------------------------------------------------------------------

/** Create identity matrix as sparse matrix */
template <typename TIdx>
DSparseMatrix<double> eye(TIdx n, TIdx procs)
{
    vector<Triplet<double, TIdx>> coefficients;
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
DSparseMatrix<double, TIdx> rand(TIdx n, TIdx m, TIdx procs, double density)
{
    vector<Triplet<double, TIdx>> coefficients;
    coefficients.reserve(n);

    double mu = 1.0 / density + 0.5;
    double sigma = 0.5 * mu;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::normal_distribution<double> gauss(mu, sigma);

    TIdx j = static_cast<int>(gauss(mt)) / 2;
    j = max(j, (TIdx)1);
    TIdx i = 0;
    while (true)
    {
        coefficients.push_back(Triplet<double, TIdx>(j, i,
                (1.0 + 10.0 * dist(mt))));

        int offset = static_cast<int>(gauss(mt));
        offset = max(offset, 1);
        j += offset;

        while (j >= n)
        {
            j = j - n;
            i++;
        }

        if (i >= n)
            break;

    }

    DSparseMatrix<double, TIdx> A(n, m);
    A.setDistributionScheme(Partitioning::cyclic, procs);

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

//-----------------------------------------------------------------------------

} // namespace Zee
