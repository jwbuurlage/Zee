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

#include "matrix/sparse.hpp"
#include "matrix/dense.hpp"
#include "parallel.hpp"

#include <vector>
using std::vector;

#include <thread>
using std::thread;

#include <iostream>
using std::cerr;
using std::endl;

namespace Zee
{

template <typename TVal, typename TIdx>
void spmv_image_cpp(DSparseMatrixImage<TVal, TIdx>& A,
        DVector<TVal>& v,
        DVector<TVal>& u)
{
    for(storage_iterator_triplets<TVal, TIdx, false> it = A.begin();
            it != A.end(); ++it) {
        u[(*it).row()] += (*it).value() * v[(*it).col()];
    }
}

template <typename TVal, typename TIdx>
void spmv_cpp(DSparseMatrix<TVal, TIdx>& A,
        DVector<TVal>& v,
        DVector<TVal>& u)
{
    vector<thread> threads;
    vector<DVector<TVal>> uis;
    for (TIdx p = 0; p < A.procs(); p++)
        uis.push_back(zeros(A.rows()));

    // start a thread with image multiplication
    TIdx proc = 0;
    for (auto& image : A.getImages())
    {
        threads.push_back(
                thread(spmv_image_cpp<TVal, TIdx>,
                    std::ref(image),
                    std::ref(v),
                    std::ref(uis[proc++])));
    }

    for(auto& t : threads)
        t.join();

    // add results (O(np), will optimize when distributing vectors)
    for (TIdx p = 0; p < A.procs(); p++)
        for (TIdx j = 0; j < A.rows(); ++j)
            u[j] += uis[p][j];
}

}
