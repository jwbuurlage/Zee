/*
File: include/partitioner.h

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include <memory>

#include "matrix/sparse.hpp"
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
            _procs = procs;
        }

    protected:
        /* The number of processors _after_ partitioning */
        int _procs = 0;

        /* The number of processors _before_ partitioning
         * The value 0 is used for an arbitrary number */
        int _procs_in = 0;
};


template <class TMatrix = DSparseMatrix<double>>
class IterativePartitioner : public Partitioner<TMatrix>
{
    public:
        IterativePartitioner() { };
        virtual ~IterativePartitioner() { };

        /** Partitioning with an IP means the initial partitioning. A call to this function is optional.
          * Default behaviour is to forward to refine. */
        virtual TMatrix& partition(
                TMatrix& A) override
        {
            return this->refine(A);
        }

        virtual TMatrix& refine(TMatrix& A) = 0;
};

enum CyclicType {
    row,
    column,
    element_wise
};

template <class TMatrix = DSparseMatrix<double>>
class CyclicPartitioner : public Partitioner<TMatrix>
{
    public:

        CyclicPartitioner()
        {
            this->_procs = 1;
        }

        CyclicPartitioner(int procs)
        {
            this->_procs = procs;
        }

        CyclicPartitioner(int procs, CyclicType type)
        {
            this->_procs = procs;
            this->_type = type;
        }

        virtual ~CyclicPartitioner() override = default;

        virtual TMatrix& partition(
                TMatrix& A) override
        {
            using TIdx = typename TMatrix::index_type;
            using TImage = typename TMatrix::image_type;

            // repartition A
            A.setDistributionScheme(Partitioning::custom, this->_procs);

            vector<unique_ptr<TImage>> new_images;
            for (TIdx i = 0; i < this->_procs; ++i)
                new_images.push_back(std::make_unique<TImage>());

            // FIXME, make this inplace for the first proc
            // FIXME: revert if/for
            auto& images = A.getMutableImages();
            TIdx cur = 0;
            for (auto& image : images)
                for(auto& trip : *image) {
                    if (_type == CyclicType::row) {
                        new_images[trip.row() % this->_procs]->pushTriplet(trip);
                    } else  if (_type == CyclicType::column) {
                        new_images[trip.col() % this->_procs]->pushTriplet(trip);
                    } else  if (_type == CyclicType::element_wise) {
                        new_images[cur++ % this->_procs]->pushTriplet(trip);
                    }

                }

            A.resetImages(new_images);

            return A;
        }

    private:
        CyclicType _type = CyclicType::row;
};

}
