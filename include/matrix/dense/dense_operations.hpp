/*
File: include/matrix/dense/dense_operations.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

// IMPLEMENTATIONS ////////////////////////////////////////////////////////////

template <typename TVal, typename TIdx>
DVector<TVal, TIdx> perform_operation(
        BinaryOperation<operation::type::addition,
        DVector<TVal, TIdx>,
        DVector<TVal, TIdx>> operation)
{
    auto& lhs = operation.getLHS();
    auto& rhs = operation.getRHS();

    ZeeAssert(lhs.size() == rhs.size());

    auto returnVector = DVector<TVal, TIdx>(lhs.size());
    for (TIdx i = 0; i < lhs.size(); ++i) {
        returnVector[i] = lhs[i] + rhs[i];
    }

    return returnVector;
}

template <typename TVal, typename TIdx>
DVector<TVal, TIdx> perform_operation(
        BinaryOperation<operation::type::subtraction,
        DVector<TVal, TIdx>,
        DVector<TVal, TIdx>> operation)
{
    auto& lhs = operation.getLHS();
    auto& rhs = operation.getRHS();

    ZeeAssert(lhs.size() == rhs.size());

    auto returnVector = DVector<TVal, TIdx>(lhs.size());
    for (TIdx i = 0; i < lhs.size(); ++i) {
        returnVector[i] = lhs[i] - rhs[i];
    }

    return returnVector;
}

template <typename TVal, typename TIdx>
DVector<TVal, TIdx> perform_operation(
        BinaryOperation<operation::type::scalar_product,
        DVector<TVal, TIdx>,
        TVal> operation)
{
    auto& lhs = operation.getLHS();
    auto& rhs = operation.getRHS();

    auto returnVector = DVector<TVal, TIdx>(lhs.size());
    for (TIdx i = 0; i < lhs.size(); ++i) {
        returnVector[i] = lhs[i] * rhs;
    }

    return returnVector;
}

template <typename TVal, typename TIdx>
DVector<TVal, TIdx> perform_operation(
        BinaryOperation<operation::type::scalar_division,
        DVector<TVal, TIdx>,
        TVal> operation)
{
    auto& lhs = operation.getLHS();
    auto& rhs = operation.getRHS();

    auto returnVector = DVector<TVal, TIdx>(lhs.size());
    for (TIdx i = 0; i < lhs.size(); ++i) {
        returnVector[i] = lhs[i] / rhs;
    }

    return returnVector;
}
