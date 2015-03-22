/*
File: include/partitioner.h

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#include "sparse_matrix.h"


namespace Zee
{
    
template <typename TVal, typename TIdx = int32_t>
class Partitioner
{
    public:
        Partitioner();
        virtual ~Partitioner();

        // Takes one sparse matrix, partitions and redistributes it.
        // TODO: Is this possible to do "in-place?"
        // TODO: What about refinement?
        virtual DSparseMatrix<TVal, TIdx> partition(
                DSparseMatrix<TVal, TIdx> A)
        {
            return -1;
        }
};

}
