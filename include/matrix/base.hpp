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
#include <memory>

#include "../operations.hpp"

namespace Zee {

/** Base class for matrices */
template <typename Derived, typename TVal, typename TIdx = int32_t>
class DMatrixBase {
    public:
        DMatrixBase(TIdx rows, TIdx cols)
            : cols_(cols),
              rows_(rows),
              procs_((TIdx)1)
        { }

        virtual ~DMatrixBase() = default;

        /** @return the total number of (possible) entries of the matrix  */
        TIdx size() const
        {
            return rows_ * cols_;
        }

        /** @return the number of rows of the matrix */
        TIdx getRows() const
        {
            return rows_;
        }

        /** @return the number of columns of the matrix */
        TIdx getCols() const
        {
            return cols_;
        }

        /** @return the number of columns of the matrix */
        TIdx getProcs() const
        {
            return procs_;
        }

        const Derived& derived() const
        {
            return *this;
        };

        #include "base_operations.hpp"

    protected:
        TIdx cols_ = 0;
        TIdx rows_ = 0;

        TIdx procs_ = 0;
};

} // namespace Zee
