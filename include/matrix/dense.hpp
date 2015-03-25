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

// FIXME: should be a specialization of a general dense matrix
template <typename TVal, typename TIdx = int32_t>
class DVector : DMatrixBase<TVal, TIdx>
{
    public:
        DVector(TIdx n) : DMatrixBase<TVal, TIdx>(n, 1)
        {}

    private:
};

DVector<double> rand(int32_t n)
{
    DVector<double> v(n);
    return v;
}

}
