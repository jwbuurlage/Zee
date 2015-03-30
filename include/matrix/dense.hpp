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

namespace Zee
{

using std::vector;

// FIXME: should be a specialization of a general dense matrix
template <typename TVal, typename TIdx = int32_t>
class DVector : public DMatrixBase<TVal, TIdx>
{
    public:
        DVector(TIdx n) : DMatrixBase<TVal, TIdx>(n, 1)
        {
            _elements.resize(n);
        }

        TVal& operator[] (int i) { return _elements[i]; }

    private:
        vector<TVal> _elements;
};

DVector<double> zeros(int32_t n)
{
    DVector<double> v(n);

    for (int32_t i = 0; i < n; ++i)
        v[i] = 0;

    return v;
}

DVector<double> rand(int32_t n)
{
    DVector<double> v(n);

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(1.0, 10.0);

    for (int32_t i = 0; i < n; ++i)
        v[i] = dist(mt);

    return v;
}

}
