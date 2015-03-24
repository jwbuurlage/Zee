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

#include <iostream>
using std::cout;
using std::endl;

#include <cstdint>

#include <vector>
using std::vector;

#include <matrix/base.hpp>
#include <matrix/dense.hpp>

namespace Zee
{

template <typename TVal, typename TIdx = int32_t>
class DSparseMatrixImage;

/** A (matrix) triplet is a tuplet (i, j, a_ij), representing an entry in a
  * matrix. It is particularly useful in the representation of sparse matrices.
  */
template <typename TVal, typename TIdx = int32_t>
class Triplet 
{
    public:
        Triplet(TIdx i, TIdx j, TIdx value)
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
enum InitialPartitioning
{
    ONE_D_CYCLIC_PARTITIONING,
    ONE_D_BLOCK_PARTITIONING,
    RANDOM_PARTITIONING
};

/** The class DSparseMatrix is a distributed matrix type inspired by 
  * Eigen's SparseMatrix. It is distributed over multiple processing units,
  */
template <typename TVal, typename TIdx = int32_t>
class DSparseMatrix : DMatrixBase<TVal, TIdx>
{
    public:
        /** Initialize an (empty) sparse (rows x cols) oatrix */
        DSparseMatrix(TIdx rows, TIdx cols) :
            DMatrixBase<TVal, TIdx>(rows, cols)
        {
            _procs = 1;
            _partitioning = ONE_D_CYCLIC_PARTITIONING;
        }

        /** Default deconstructor */
        ~DSparseMatrix() { }

        /** Move constructor
          * FIXME: COPY RESOURCES
          */
        DSparseMatrix(DSparseMatrix&& o) = default;

        /** Sets the distribution scheme for this matrix */
        void setDistributionScheme(InitialPartitioning partitioning,
                TIdx procs)
        {
            _partitioning = partitioning;
            _procs = procs;
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
            // TODO: implement
            return -1;
        }

        /** Construct a matrix from a set of triplets */
        template<typename TInputIterator>
        void setFromTriplets(
            const TInputIterator& begin,
            const TInputIterator& end)
        {
            if(!_subs.empty()) {
                _subs.clear();
                _nz = 0;
            }

            //TODO: implement
            for (TInputIterator it = begin; it != end; it++)
            {
                cout << "(" << (*it).col() <<
                    ", " << (*it).row() <<
                    ", " << (*it).value() << ")" << endl;

                _nz++;
            }

            return;
        }

    private:
        TIdx _nz;

        // TODO: how do we reference other processors that store (part) of spm
        TIdx _procs;
        InitialPartitioning _partitioning;
        vector<DSparseMatrixImage<TVal, TIdx>> _subs;
};

/** Storage type for sparse matrix (image). */
enum StorageType
{
    COMPRESSED_ROW,
    COMPRESSED_COLUMN
};

// created and owned by a processor. It is a submatrix, which holds the actual
// data, the 'global' DSparseMatrix can be seen as the sum of these images.
template <typename TVal, typename TIdx>
class DSparseMatrixImage
{
    public:
        DSparseMatrixImage() { }
        ~DSparseMatrixImage() { }

    private:
        // Compressed Row Storage or Compressed Column Storage
};

//-----------------------------------------------------------------------------
// Convenience functions (MATLAB syntax)
//-----------------------------------------------------------------------------

// Create identity matrix as sparse matrix
DSparseMatrix<double> eye(int32_t n)
{
    // construct coefficients
    vector<Triplet<double>> coefficients;
    coefficients.reserve(n);

    for (int i = 0; i < n; ++i)
        coefficients.push_back(Triplet<double>(i, i, 1.0));

    DSparseMatrix<double> A(n, n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

DSparseMatrix<double> rand(int32_t n, int32_t m, double fill_in)
{
    // TODO: fill with random stuff, use storage

    vector<Triplet<double>> coefficients;
    coefficients.reserve(n);

    for (int i = 0; i < n; ++i)
        coefficients.push_back(Triplet<double>(i, i, 1.0));

    DSparseMatrix<double> A(n, m);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

//-----------------------------------------------------------------------------

}
