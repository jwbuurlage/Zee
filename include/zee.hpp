/*
File: include/zee.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include "default_types.hpp"

#include "matrix/base/base.hpp"
#include "matrix/sparse/sparse.hpp"
#include "matrix/dense/dense.hpp"

#include "operations/operation_types.hpp"
#include "operations/operations.hpp"

#include "partitioners/partitioner.hpp"
#include "partitioners/kernighan_lin.hpp"
#include "partitioners/medium_grain.hpp"
#include "partitioners/multi_level.hpp"
#include "partitioners/vector_partitioner.hpp"

#include "matrix_market.hpp"
#include "matrix_toolbox.hpp"

#include "benchmarking.hpp"
#include "logging.hpp"
#include "plotter.hpp"
#include "common.hpp"
