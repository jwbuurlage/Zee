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
        using TIdx = typename TMatrix::index_type;

        /* Use the vector distribution to compute local indices and
         * propagate this to storage */
        // MAKE THIS GENERAL, do not assume dist(u) = dist(v)
        // 1. image should receive a list of its vector indices I_s { j | P(v_j) = s }
        // 2. image then constructs a map i = {1, 2, 3, ..} -> j -> I_s
        // 3. storage updates using the inverse of this map
        TIdx p = A_.getProcs();

        // 1. construct list of vector indices
        std::vector<std::vector<TIdx>> localIndicesV(p);
        auto& ownersV = v_.getOwners();
        TIdx i = 0;
        for (auto owner : ownersV) {
            localIndicesV[owner].push_back(i++);
        }

        std::vector<std::vector<TIdx>> localIndicesU(p);
        auto& ownersU = u_.getOwners();
        i = 0;
        for (auto owner : ownersU) {
            localIndicesU[owner].push_back(i++);
        }

        // 2 + 3. inform images
        TIdx s = 0;
        for (auto& image : A_.getMutableImages()) {
            image->setLocalIndices(std::move(localIndicesV[s]),
                                   std::move(localIndicesU[s]));
            s++;
        }
        localIndicesV.clear();
        localIndicesU.clear();

        // set owners for images
        for (auto& image : A_.getMutableImages()) {
            for (TIdx idx = image->getNumLocalV();
                 idx < image->getLocalIndicesV().size(); ++idx) {
                image->getRemoteOwnersV().push_back(
                    v_.getOwners()[image->getLocalIndicesV()[idx]]);
            }
            for (TIdx idx = image->getNumLocalU();
                 idx < image->getLocalIndicesU().size(); ++idx) {
                image->getRemoteOwnersU().push_back(
                    u_.getOwners()[image->getLocalIndicesU()[idx]]);
            }

            image->localizeStorage();
        }
    }

  protected:
    TMatrix& A_;
    TVector& v_;
    TVector& u_;
};

} // namespace Zee
