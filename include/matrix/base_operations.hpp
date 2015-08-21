/*
File: include/matrix/base_operations.hpp

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

// OPERATORS //////////////////////////////////////////////////////////////////

BinaryOperation<operation::type::scalar_product,
    Derived, TVal>
    operator*(const TVal& rhs) const
{
    return BinaryOperation<operation::type::scalar_product,
           Derived, TVal>(*this, rhs);
}

/** Scalar multiplication */
BinaryOperation<operation::type::scalar_division,
    Derived, TVal>
    operator/(const TVal& rhs) const
{
    return BinaryOperation<operation::type::scalar_division,
           Derived, TVal>(*this, rhs);
}

// TODO: template OtherDerived
BinaryOperation<operation::type::subtract,
    Derived, Derived>
    operator- (const Derived& rhs) const
{
    return BinaryOperation<operation::type::subtract,
           Derived, Derived>(derived(), rhs);
}

template <operation::type TOp, typename TLHS, typename TRHS>
BinaryOperation<operation::type::subtract,
    Derived, Derived>
 operator- (const BinaryOperation<TOp, TLHS, TRHS>& op) const
{
    Derived v(this->size());
    v = op;
    ZeeLogVar(v.size());
    auto subOp = *this - v;
    return subOp;
}

BinaryOperation<operation::type::add,
    Derived, Derived>
    operator+ (const Derived& rhs) const
{
    return BinaryOperation<operation::type::add,
           Derived, Derived>(*this, rhs);
}

template <operation::type TOp, typename TLHS, typename TRHS>
BinaryOperation<operation::type::add,
    Derived, Derived>
 operator+ (const BinaryOperation<TOp, TLHS, TRHS>& op) const
{
    Derived v(this->size());
    v = op;
    ZeeLogVar(v.size());
    auto addOp = *this + v;
    return addOp;
}

template <operation::type TOp, typename TLHS, typename TRHS>
void operator-= (const BinaryOperation<TOp, TLHS, TRHS>& op)
{
    Derived v(this->size());
    v = op;
    auto subOp = *this - v;
    *this = subOp;
}

void operator-= (const Derived& rhs)
{
    auto op = *this - rhs;
    *this = op;
}
