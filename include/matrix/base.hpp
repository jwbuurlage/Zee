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

#include <cstdint>
#include <vector>
using std::vector;


namespace Zee
{

/** Base class for matrices */
template <typename TVal, typename TIdx = int32_t>
class DMatrixBase
{
    public:
        DMatrixBase(TIdx rows, TIdx cols)
        {
            _rows = rows;
            _cols = cols;
        }

        virtual ~DMatrixBase() { };

        /** @return the total number of (possible) entries of the matrix  */
        inline TIdx size() const { return _rows * _cols; }

        /** @return the number of rows of the matrix */
        inline TIdx rows() const { return _rows; }

        /** @return the number of columns of the matrix */
        inline TIdx cols() const { return _cols; }

    protected:
        TIdx _rows;
        TIdx _cols;
};

}
