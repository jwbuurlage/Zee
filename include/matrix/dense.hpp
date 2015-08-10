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

#include <unpain_base.hpp>

#include "base.hpp"
#include "sparse.hpp"
#include "../operations.hpp"

namespace Zee {

using std::vector;

// FIXME: should be a specialization of a general dense matrix
// FIXME: saved as pairs? or just owners distributed cyclically
template <typename TVal, typename TIdx = int32_t>
class DVector : public DMatrixBase<TVal, TIdx>
{
    public:
        explicit DVector(std::shared_ptr<UnpainBase::Center<TIdx>> center, TIdx n, TVal defaultValue = (TVal)0)
            : DMatrixBase<TVal, TIdx>(center, n, 1)
        {
            _elements.resize(n);
            std::fill(_elements.begin(), _elements.end(), defaultValue);
        }

        TVal& operator[] (int i)
        {
            return _elements[i];
        }

        const TVal& operator[] (int i) const
        {
            return _elements[i];
        }

        void operator= (BinaryOperation<operation::type::product,
            DSparseMatrix<TVal, TIdx>,
            DVector<TVal, TIdx>> op)
        {
            auto& A = op.getLHS();
            auto& v = op.getRHS();
            auto& u = *this;

            // FIXME: parallelize
            for (auto& submatrixPtr : A.getImages()) {
                for (auto& triplet : *submatrixPtr) {
                    u[triplet.row()] = triplet.value() * v[triplet.col()];
                }
            }
        }

    private:
        vector<TVal> _elements;
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
