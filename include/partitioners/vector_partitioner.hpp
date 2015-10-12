/*
File: include/partitioner.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

// FIXME ?? dont want this just predeclarations
#include "matrix/sparse/sparse.hpp"
#include "matrix/dense/dense.hpp"

namespace Zee {

template <class TMatrix = DSparseMatrix<>, class TVector = DVector<>>
class VectorPartitioner {
  public:
    VectorPartitioner(TMatrix& A, TVector& v, TVector& u)
        : A_(A), v_(v), u_(u) {}

    virtual void partition() {}

    virtual void localizeMatrix() {
        /* Use the vector distribution to compute local indices and
         * propagate this to storage */
        // 1. image should receive a list of its vector indices I_s { j | P(v_j) = s }
        // 2. image then constructs a map i = {1, 2, 3, ..} -> j -> I_s
        // 3. storage updates using the inverse of this map
        // 4. MAKE THIS GENERAL, do not assume dist(u) = dist(v)
        ZeeLogDebug << "Localizing" << endLog;
    }

  protected:
    TMatrix& A_;
    TVector& v_;
    TVector& u_;
};

} // namespace Zee
