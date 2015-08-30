/*
File: include/matrix/base/base.hpp

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

#include "../../default_types.hpp"
#include "../../operations/operation_types.hpp"

namespace Zee {

template <operation::type op, typename TLHS, typename TRHS>
class BinaryOperation;

/** Base class for matrices */
template <typename Derived,
         typename TVal = default_scalar_type,
         typename TIdx = default_index_type>
class DMatrixBase {
    public:
        DMatrixBase(TIdx rows, TIdx cols)
            : cols_(cols),
              rows_(rows),
              procs_(1)
        { }

        DMatrixBase()
            : DMatrixBase(0, 0)
        { }

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

        Derived& derived()
        {
            return *(Derived*)this;
        };

        const Derived& derived() const
        {
            return *(Derived*)this;
        };

        #include "base_operations.hpp"

    protected:
        TIdx cols_ = 0;
        TIdx rows_ = 0;

        TIdx procs_ = 0;
};

#include "base_operations_global.hpp"
} // namespace Zee
