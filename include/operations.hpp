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

#include "matrix/sparse.hpp"
#include "matrix/dense.hpp"
#include "parallel.hpp"

namespace Zee
{

template<ParallelProvider P>
void SpMV<P>(DSparseMatrix& A, DVector& v, DVector& u);

void SpMV<P_CPP>(DSparseMatrix& A,DVector& v, DVector& u)
{
    return;
}

}
