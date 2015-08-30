/*
File: include/matrix/base/base_operations.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

// FIXME: IIRC boost has support to wrap this abstractly with ops which would
// minimize the LOC needed for these operations
//
// Here we mimic the behaviour of Eigen
// http://eigen.tuxfamily.org/dox/TopicInsideEigenExample.html

template <operation::type S, typename T, typename U>
void operator= (const BinaryOperation<S, T, U>& op) {
    derived() = std::move(perform_operation(op));
}

template <operation::type S, typename T, typename U>
void operator-= (const BinaryOperation<S, T, U>& op) {
    auto w = perform_operation(op);
    derived() = derived() - w;
}

template <operation::type S, typename T, typename U>
void operator+= (const BinaryOperation<S, T, U>& op) {
    auto w = perform_operation(op);
    derived() = derived() + w;
}

template <operation::type S, typename T, typename U>
void operator*= (const BinaryOperation<S, T, U>& op) {
    auto w = perform_operation(op);
    derived() = derived() * w;
}

template <operation::type S, typename T, typename U>
void operator/= (const BinaryOperation<S, T, U>& op) {
    auto w = perform_operation(op);
    derived() = derived() / w;
}

template <typename OtherDerived, typename S, typename T>
BinaryOperation<operation::type::addition, Derived, OtherDerived>
operator+ (const DMatrixBase<OtherDerived, S, T>& other) const
{
    return BinaryOperation<operation::type::addition,
           Derived,
           OtherDerived>(derived(), other.derived());
}

template <typename OtherDerived, typename S, typename T>
BinaryOperation<operation::type::subtraction, Derived, OtherDerived>
operator- (const DMatrixBase<OtherDerived, S, T>& other) const
{
    return BinaryOperation<operation::type::subtraction,
           Derived,
           OtherDerived>(derived(), other.derived());
}

template <typename OtherDerived, typename S, typename T>
BinaryOperation<operation::type::product, Derived, OtherDerived>
operator* (const DMatrixBase<OtherDerived, S, T>& other) const
{
    return BinaryOperation<operation::type::product,
           Derived,
           OtherDerived>(derived(), other.derived());
}

// Maybe ideal candidate for Concepts (See TS or C++17)
// problematic when doing alpha * (x + y) where x, y vectors?
// we need to match scalar types
BinaryOperation<operation::type::scalar_product, Derived, TVal>
operator* (const TVal& alpha) const
{
    return BinaryOperation<operation::type::scalar_product,
           Derived,
           TVal>(derived(), alpha);
}

BinaryOperation<operation::type::scalar_division, Derived, TVal>
operator/ (const TVal& alpha) const
{
    return BinaryOperation<operation::type::scalar_division,
           Derived,
           TVal>(derived(), alpha);
}
