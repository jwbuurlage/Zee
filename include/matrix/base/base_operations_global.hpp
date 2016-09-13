/*
File: include/matrix/base/base_operations_global.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

// Allow expression of type lambda * V by forwarding to V * lambda
template <typename Derived, typename TVal, typename TIdx>
BinaryOperation<operation::type::scalar_product, Derived, TVal> operator*(
    const TVal& alpha, const DMatrixBase<Derived, TVal, TIdx>& base) {
    return BinaryOperation<operation::type::scalar_product, Derived, TVal>(
        base.derived(), alpha);
}
