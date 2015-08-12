/*
File: include/matrix/dense_operations.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

// OPERATORS //////////////////////////////////////////////////////////////////

/** Scalar multiplication */
BinaryOperation<operation::type::scalar_product,
    DVector<TVal, TIdx>, TVal>
    operator*(const TVal& rhs) const
{
    return BinaryOperation<operation::type::scalar_product,
           DVector<TVal, TIdx>, TVal>(*this, rhs);
}

/** Scalar multiplication */
BinaryOperation<operation::type::scalar_division,
    DVector<TVal, TIdx>, TVal>
    operator/(const TVal& rhs) const
{
    return BinaryOperation<operation::type::scalar_division,
           DVector<TVal, TIdx>, TVal>(*this, rhs);
}

BinaryOperation<operation::type::subtract,
    DVector<TVal, TIdx>, DVector<TVal, TIdx>>
    operator- (const DVector<TVal, TIdx>& rhs)
{
    return BinaryOperation<operation::type::subtract,
           DVector<TVal, TIdx>, DVector<TVal, TIdx>>(*this, rhs);
}

BinaryOperation<operation::type::add,
    DVector<TVal, TIdx>, DVector<TVal, TIdx>>
    operator+ (const DVector<TVal, TIdx>& rhs)
{
    return BinaryOperation<operation::type::add,
           DVector<TVal, TIdx>, DVector<TVal, TIdx>>(*this, rhs);
}

template <operation::type TOp, typename TLHS, typename TRHS>
void operator-= (const BinaryOperation<TOp, TLHS, TRHS>& op)
{
    DVector<TVal, TIdx> v(this->getCenter(), this->size());
    v = op;
    auto subOp = *this - v;
    *this = subOp;
}

void operator-= (const DVector<TVal, TIdx>& rhs)
{
    auto op = *this - rhs;
    *this = op;
}

// IMPLEMENTATIONS ////////////////////////////////////////////////////////////

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

// We cannot rely soleley on scalar_product because we pass rhs by reference
// here we let the reciprocal live long enough for the assignment to take place
void operator= (BinaryOperation<operation::type::scalar_division,
    DVector<TVal, TIdx>, TVal> op)
{
    const auto& v = op.getLHS();
    const auto& lambda = op.getRHS();
    TVal reciprocal = (TVal)1 / lambda;
    auto scalarProductOp = v * reciprocal;
    *this = scalarProductOp;
}

void operator= (BinaryOperation<operation::type::scalar_product,
    DVector<TVal, TIdx>, TVal> op)
{
    const auto& v = op.getLHS();
    const auto& lambda = op.getRHS();

    if (v.size() != size()) {
        ZeeLogError << "Can not assign DVector to expression of type"
            " (alpha * DVector) with a different dimension." << endLog;
        return;
    }

    for (TIdx i = 0; i < this->size(); ++i) {
        elements_[i] = lambda * v[i];
    }
}

void operator= (BinaryOperation<operation::type::subtract,
    DVector<TVal, TIdx>, DVector<TVal, TIdx>> op)
{

    const auto& lhs = op.getLHS();
    const auto& rhs = op.getRHS();

    if (rhs.size() != lhs.size()) {
        ZeeLogError << "Can subtract vectors of different sizes" << endLog;
        return;
    } else if (lhs.size() != size()) {
        ZeeLogError << "Cannot assign DVector to expression of type"
            " (v - w) with a different dimension." << endLog;
        return;
    }

    for (TIdx i = 0; i < this->size(); ++i) {
        elements_[i] = lhs[i] - rhs[i];
    }
}
