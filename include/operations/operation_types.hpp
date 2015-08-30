/*
File: include/operations/operation_types.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

namespace Zee {

namespace operation
{
    enum class type
    {
        product,
        addition,
        subtraction,
        scalar_product,
        scalar_division
    };
} // namespace operation

} // namespace Zee
