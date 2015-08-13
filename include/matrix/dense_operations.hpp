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
    const auto p = A.getProcs();

    std::mutex writeMutex;
    Barrier<TIdx> barrier(p);

    // FIXME: parallelize using unpain
    A.compute([&] (std::shared_ptr<TImage> submatrixPtr,
                TIdx s) {
        std::map<TIdx, TVal> u_s;
        for (const auto& triplet : *submatrixPtr) {
             u_s[triplet.row()] += triplet.value() * v[triplet.col()];
        }

        // if RHS = this then we can only reset here,
        // so we use a barrier sync
        barrier.sync();

        for (TIdx i = s; i < this->size(); i += p) {
            u[i] = 0;
        }

        barrier.sync();

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

    ZeeAssert(v.size() == size());

    for (TIdx i = 0; i < this->size(); ++i) {
        elements_[i] = lambda * v[i];
    }
}

void operator= (BinaryOperation<operation::type::subtract,
    DVector<TVal, TIdx>, DVector<TVal, TIdx>> op)
{
    const auto& lhs = op.getLHS();
    const auto& rhs = op.getRHS();

    ZeeAssert(rhs.size() == lhs.size());
    ZeeAssert(lhs.size() == size());

    for (TIdx i = 0; i < this->size(); ++i) {
        elements_[i] = lhs[i] - rhs[i];
    }
}
