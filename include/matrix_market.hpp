/*
File: include/matrix_market.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

// # TODO
// [ ] Support the extended format

#pragma once

#include <fstream>
#include <string>
#include <sstream>

#include "logging.hpp"

namespace Zee {

namespace matrix_market
{
    enum info : int {
        coordinate     = 1 << 0,
        real           = 1 << 1,
        pattern        = 1 << 2,
        general        = 1 << 3,
        symmetric      = 1 << 4,
        skew_symmetric = 1 << 5
    };
} //namespace matrix_market

} // namespace Zee
