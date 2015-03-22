/*
File: include/sparse_matrix.h

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#include <cstdint>
#include <vector>
using std::vector;


namespace Zee
{

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

        ~Triplet();

        /** @return the value for this triplet */
        inline TVal value() const
        {
            return _value;
        };

        /** @return the row position inside the matrix */
        inline TIdx row() const
        {
            return _i;
        };

        /** @return the column position inside the matrix */
        inline TIdx col() const
        {
            return _j;
        };

    private:
        TIdx _i;
        TIdx _j;
        TVal _value;
};

/** The class DSparseMatrix is a distributed matrix type inspired by 
  * Eigen's SparseMatrix. It is distributed over multiple processing units,
  */
template <typename TVal, typename TIdx = int32_t>
class DSparseMatrix
{
    public:
        DSparseMatrix(TIdx rows, TIdx cols)
        {
            _rows = rows;
            _cols = cols;
        };

        ~DSparseMatrix();

        //---------------------------------------------------------------------
        // Information on size of matrix
        //---------------------------------------------------------------------
        /** @return the number of rows of the matrix */
        inline TIdx rows() const
        {
            return _rows;
        };

        /** @return the number of columns of the matrix */
        inline TIdx cols() const
        {
            return _cols;
        };

        /** @return the total size of the matrix (rows * cols) */
        inline TIdx size() const
        {
            return _rows * _cols;
        };

        //---------------------------------------------------------------------
        // Information on sparsity of matrix
        //---------------------------------------------------------------------
        TIdx nonZeros() const
        {
            // TODO: implement
            return -1;
        }

        void setFromTriplets(const iterator<Triplet> begin,
            const iterator<Triplet> end)
        {

        }

    private:
        TIdx _rows;
        TIdx _cols;

        vector< Triplet<TVal, TIdx> > _elements;
        // TODO: how do we reference other processors that store (part) of spm
        //vector< > _siblings;
};

// Create identity matrix
DSparseMatrix<double>* eye(int32_t n)
{
    // construct coefficients
    vector< Triplet<double> > coefficients;
    coefficients.reserve(n);

    for (int i = 0; i < n; ++i) {
        coefficients.push_back(Triplet(i, i, 1.0));
    }

    DSparseMatrix A(n, n);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

}
