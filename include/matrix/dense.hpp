/*
File: include/matrix/dense.hpp

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
#include <random>
#include <memory>
#include <atomic>
#include <mutex>
#include <map>

#include "base.hpp"
#include "common.hpp"
#include "sparse.hpp"
#include "../operations.hpp"

namespace Zee {

//-----------------------------------------------------------------------------
// VECTOR

// FIXME: should be a specialization of a general dense matrix
// FIXME: saved as pairs? or just owners distributed cyclically
template <typename TVal, typename TIdx = int32_t>
class DVector : public DMatrixBase<DVector<TVal, TIdx>, TVal, TIdx>
{
    using Base = DMatrixBase<DVector<TVal, TIdx>, TVal, TIdx>;

    public:
        explicit DVector(TIdx n, TVal defaultValue = (TVal)0)
            : Base(n, 1)
        {
            elements_.resize(n);
            std::fill(elements_.begin(), elements_.end(), defaultValue);
        }

        /** Return the size of the matrix (i.e. the number of rows) */
        TIdx size() const {
            return this->rows_;
        }

        TVal& operator[] (TIdx i)
        {
            return elements_[i];
        }

        const TVal& operator[] (TIdx i) const
        {
            return elements_[i];
        }

        /** The dot product */
        // FIXME parallellize
        TVal dot(const DVector<TVal, TIdx>& rhs) const {
            ZeeAssert(rhs.size() == size());

            const auto& lhs = *this;

            TVal sum = 0;
            for (TIdx i = 0; i < this->size(); ++i) {
                sum += lhs[i] * rhs[i];
            }

            return sum;
        }

        /** The (Euclidian) vector norm */
        TVal norm() const {
            TVal sum = 0;
            for (auto& element : elements_) {
                sum += element * element;
            }
            return sqrt(sum);
        }

        void reset() {
            for (auto& elem : elements_) {
                elem = 0;
            }
        }

        // Operator overloads and algorithm implementations
        #include "dense_operations.hpp"

    private:
        std::vector<TVal> elements_;
};

// VECTOR LOGGING

template <typename TVal, typename TIdx>
Logger& operator <<(Logger& lhs, const DVector<TVal, TIdx>& rhs) {
    auto sep = "";
    lhs << "[";
    for (TIdx i = 0; i < rhs.size(); ++i) {
        lhs << std::fixed << std::setprecision(2) << sep << rhs[i];
        sep = ", ";
    }
    lhs << "]";

    return lhs;
}

// Allow expression of type lambda * V by forwarding to V * lambda
template <typename TVal, typename TIdx>
BinaryOperation<operation::type::scalar_product,
    DVector<TVal, TIdx>, TVal>
    operator*(const TVal& lhs, const DVector<TVal, TIdx>& rhs)
{
    return BinaryOperation<operation::type::scalar_product,
           DVector<TVal, TIdx>, TVal>(rhs, lhs);
}

//-----------------------------------------------------------------------------
// MATRIX

template <typename TVal, typename TIdx = int32_t>
class DMatrix : public DMatrixBase<DMatrix<TVal, TIdx>, TVal, TIdx>
{
    public:
       // Operator overloads and algorithm implementations
        //#include "dense_matrix_operations.hpp"

    private:
        std::vector<std::vector<TVal>> elements_;
};

} // namespace Zee
