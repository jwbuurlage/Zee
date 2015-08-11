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

        // FIXME we can use #include to split this class over multiple files
        // which would make more sense
        void operator= (BinaryOperation<operation::type::product,
            DSparseMatrix<TVal, TIdx>,
            DVector<TVal, TIdx>> op)
        {
            using TImage = typename DSparseMatrix<TVal, TIdx>::image_type;

            const auto& A = op.getLHS();
            const auto& v = op.getRHS();
            auto& u = *this;

            std::mutex writeMutex;

            // FIXME: parallelize using unpain
            A.compute([&v, &u, &writeMutex] (std::shared_ptr<TImage> submatrixPtr) {
                std::map<TIdx, TVal> u_s;
                for (const auto& triplet : *submatrixPtr) {
                    u_s[triplet.row()] = u_s[triplet.row()] + triplet.value() * v[triplet.col()];
                }

                for (const auto& pKeyValue : u_s) {
                    std::lock_guard<std::mutex> lock(writeMutex);
                    u[pKeyValue.first] += pKeyValue.second;
                }
            });
        }

        /** The dot product */
        // FIXME parallellize
        TVal dot(const DVector<TVal, TIdx>& rhs) const {
            if (rhs.size() != size()) {
                ZeeLogError << "Can not compute dotproduct between vectors of"
                    "different sizes" << endLog;
                return (TVal)0;
            }

            TVal result = (TVal)0;
        }

        /** The (Euclidian) vector norm */
        TVal norm() const {
            return (TVal)0;
        }

    private:
        std::vector<TVal> elements_;
};

//-----------------------------------------------------------------------------

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
