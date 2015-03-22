/*
File: include/sparse_matrix.h

This file is part of the Zee partitioning framework

Copyright (C) 2014 Jan-Willem Buurlage
Support e-mail: <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
and the GNU Lesser General Public License along with this program,
see the file LICENSE. If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdint>
#include <vector>
using std::vector;


namespace Zee {

/* A (matrix) triplet is a tuplet (i, j, a_ij), representing an entry in a
 * matrix. It is particularly useful in the representation of sparse matrices.
 */
template <typename TVal, typename TIdx = int32_t>
class Triplet 
{
    public:
        Triplet();

    private:
        TIdx _i;
        TIdx _j;
        TVal _val;
}

/* The class DSparseMatrix is a distributed matrix type inspired by 
 * Eigen's SparseMatrix. It is distributed over multiple processing units,
 */
template <typename TVal, typename TIdx = int32_t>
class DSparseMatrix
{
    public:
        DSparseMatrix();
        ~DSparseMatrix();

    private:
        TIdx _rows;
        TIdx _cols;

        vector< Triplet<TVal, TIdx> > _elements;
        // TODO: how do we reference other processors that store (part) of spm
        //vector< > _siblings;
}

}
