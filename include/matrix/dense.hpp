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
#include <random>
#include <memory>
#include <atomic>
#include <mutex>
#include <map>

#include <unpain_base.hpp>

#include "base.hpp"
#include "common.hpp"
#include "sparse.hpp"
#include "../operations.hpp"

namespace Zee {

// FIXME: should be a specialization of a general dense matrix
// FIXME: saved as pairs? or just owners distributed cyclically
template <typename TVal, typename TIdx = int32_t>
class DVector : public DMatrixBase<TVal, TIdx>
{
    public:
        explicit DVector(std::shared_ptr<UnpainBase::Center<TIdx>> center, TIdx n, TVal defaultValue = (TVal)0)
            : DMatrixBase<TVal, TIdx>(center, n, 1)
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
            if (rhs.size() != size()) {
                ZeeLogError << "Can not compute dotproduct between vectors of"
                    "different sizes" << endLog;
                return 0;
            }

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

//-----------------------------------------------------------------------------
// LOGGING
//     Vector2D Vector2D::operator+(const Vector2D& right)const {...}

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

// FIXME: these should be constructors
//DVector<double> zeros(int32_t n)
//{
//    DVector<double> v(n);
//
//    for (int32_t i = 0; i < n; ++i)
//        v[i] = 0;
//
//    return v;
//}
//
//DVector<double> rand(int32_t n)
//{
//    DVector<double> v(n);
//
//    std::random_device rd;
//    std::mt19937 mt(rd());
//    std::uniform_real_distribution<double> dist(1.0, 10.0);
//
//    for (int32_t i = 0; i < n; ++i)
//        v[i] = dist(mt);
//
//    return v;
//}

} // namespace Zee
