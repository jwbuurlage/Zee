/*
File: include/operations/operations.hpp

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

#include "operation_types.hpp"
#include "../matrix/base/base.hpp"

namespace Zee {

template <operation::type opType, typename LHS, typename RHS>
class BinaryOperation
    : public DMatrixBase<BinaryOperation<opType, LHS, RHS>>
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

template <operation::type optype_1,
          operation::type optype_2,
          typename S, typename T, typename U>
auto perform_operation(BinaryOperation<optype_1, S,
        BinaryOperation<optype_2, T, U>> op) {
    return perform_operation(
            BinaryOperation<optype_1, S,
                decltype(perform_operation(op.getRHS()))>(
                op.getLHS(),
                perform_operation(op.getRHS()))
        );
};

template <operation::type optype_1,
          operation::type optype_2,
          typename S, typename T, typename U>
auto perform_operation(BinaryOperation<optype_1, BinaryOperation<optype_2, S, T>, U> op) {
    return perform_operation(
            BinaryOperation<optype_1,
                decltype(perform_operation(op.getLHS())), U>(
                perform_operation(op.getLHS()),
                op.getRHS())
        );
};

template <operation::type optype_1,
          operation::type optype_2,
          operation::type optype_3,
          typename S, typename T, typename U, typename V>
auto perform_operation(BinaryOperation<optype_1,
        BinaryOperation<optype_2, S, T>,
        BinaryOperation<optype_3, U, V>> op) {
    return perform_operation(
            BinaryOperation<optype_1,
                decltype(perform_operation(op.getLHS())),
                decltype(perform_operation(op.getRHS()))>(
                    perform_operation(op.getLHS()),
                    perform_operation(op.getRHS()))
            );
};

} // namespace Zee
