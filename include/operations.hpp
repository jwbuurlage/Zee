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

#include <vector>
#include <thread>
#include <iostream>

namespace Zee {

namespace operation
{
    enum type
    {
        product,
        sum
    };
} // namespace operation

template <operation::type opType, typename LHS, typename RHS>
class BinaryOperation
{
    public:
        BinaryOperation(const LHS& lhs, const RHS& rhs)
            : lhs_(lhs), rhs_(rhs)
        {
        }

        const LHS& getLHS() const {
            return lhs_;
        }

        const RHS& getRHS() const {
            return rhs_;
        }

    private:
        const LHS& lhs_;
        const RHS& rhs_;
};

} // namespace Zee
