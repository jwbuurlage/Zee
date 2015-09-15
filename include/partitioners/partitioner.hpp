/*
File: include/partitioner.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include <memory>

#include "matrix/sparse/sparse.hpp"
#include "common.hpp"

namespace Zee {

using std::unique_ptr;
using std::vector;

/**
 * This class (re)partitions a sparse matrix.
 * FIXME: CRT for factory?
 */
template <class TMatrix = DSparseMatrix<double>>
class Partitioner
{
    public:
        Partitioner() { }

        virtual ~Partitioner() { };

        // Takes one sparse matrix, partitions and redistributes it.
        // TODO: Is this possible to do "in-place?"
        virtual TMatrix& partition(
                TMatrix& A) = 0;

        virtual void initialize(TMatrix& A) { }

        /** Set the number of processors
          * @param procs Number of processors _after_ partitioning.
          * */
        void setProcs(int procs) {
            procs_ = procs;
        }

    protected:
        /* The number of processors _after_ partitioning */
        int procs_ = 0;

        /* The number of processors _before_ partitioning
         * The value 0 is used for an arbitrary number */
        int procs_in_ = 0;
};


template <class TMatrix = DSparseMatrix<double>>
class IterativePartitioner : public Partitioner<TMatrix>
{
    public:
        // FIXME should probably hold a shared_ptr to A, dangling reference danger
        IterativePartitioner(TMatrix& A) : A_(A) {};
        virtual ~IterativePartitioner() { };

        /** Partitioning with an IP means the initial partitioning. A call to this function is optional.
          * Default behaviour is to forward to refine. */
        virtual void partition(
                TMatrix& A) override
        {
            return this->refine();
        }

        virtual void refine() = 0;

    protected:
        TMatrix& A_;
};

} // namespace Zee
