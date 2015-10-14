/*
File: include/matrix/dense/dense_matrix_operations.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

template <typename TVal, typename TIdx>
DVector<TVal, TIdx> perform_operation(
    BinaryOperation<operation::type::product, DSparseMatrix<TVal, TIdx>,
                    DVector<TVal, TIdx>> op) {
    using TImage = typename DSparseMatrix<TVal, TIdx>::image_type;

    const auto& A = op.getLHS();
    const auto& v = op.getRHS();
    DVector<TVal, TIdx> u{A.getRows(), 0.0};
    const auto p = A.getProcs();

    ZeeAssert(A.getCols() == v.size());

    std::mutex writeMutex;
    Barrier<TIdx> barrier(p);

    A.compute([&](std::shared_ptr<TImage> submatrixPtr, TIdx s) {
        std::vector<TIdx> localU(submatrixPtr->getLocalIndicesU().size());
        auto& localIndicesV = submatrixPtr->getLocalIndicesV();
        for (const auto& triplet : *submatrixPtr) {
            // FIXME V should be distributed too
            localU[triplet.row()] +=
                triplet.value() * v[localIndicesV[triplet.col()]];
        }

        // if RHS is a reference to *this then we can only reset here,
        // so we use a barrier sync
        barrier.sync();

        for (TIdx i = s; i < u.size(); i += p) {
            u[i] = 0;
        }

        barrier.sync();

        TIdx i = 0;
        auto& localIndicesU = submatrixPtr->getLocalIndicesU();
        for (const auto& val : localU) {
            std::lock_guard<std::mutex> lock(writeMutex);
            u[localIndicesU[i++]] += val;
        }
    });

    return u;
}

template <typename TVal, typename TIdx>
DMatrix<TVal, TIdx>
perform_operation(BinaryOperation<operation::type::product, DMatrix<TVal, TIdx>,
                                  DMatrix<TVal, TIdx>> op) {
    const auto& A = op.getLHS();
    const auto& B = op.getRHS();
    DMatrix<TVal, TIdx> C(A.getRows(), B.getCols());

    ZeeAssert(A.getCols() == B.getRows());

    for (TIdx i = 0; i < C.getRows(); ++i) {
        for (TIdx j = 0; j < C.getRows(); ++j) {
            for (TIdx k = 0; k < A.getCols(); ++k) {
                C.at(i, j) += A.at(i, k) * B.at(k, j);
            }
        }
    }

    return C;
}
