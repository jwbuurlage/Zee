/*
File: include/partitioner.h

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#include <memory>

#include "matrix/sparse.hpp"

namespace Zee {

using std::unique_ptr;
using std::vector;

template<class TPart>
class Factory
{
    public:
        Factory() { };
        ~Factory() { };

        unique_ptr<TPart> make()
        {
            return unique_ptr<TPart>(new TPart());
        }
};

/** 
 * This class (re)partitions a sparse matrix. 
 * FIXME: CRT for factory?
 */
template <class TMatrix = DSparseMatrix<double>>
class Partitioner
{
    public:
        Partitioner() { };
        virtual ~Partitioner() { };

        // Takes one sparse matrix, partitions and redistributes it.
        // TODO: Is this possible to do "in-place?"
        // TODO: What about refinement?
        virtual TMatrix& partition(
                TMatrix& A) = 0;

        /** Set the number of processors
          * @param procs Number of processors _after_ partitioning.
          * */
        void setProcs(int procs)
        {
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
class IterativePartitioner : Partitioner<TMatrix>
{
    public:
        IterativePartitioner() { };
        virtual ~IterativePartitioner() { };

        /** Partitioning with an IP is simply refining. This function
         * forwards to refine. */
        virtual TMatrix& partition(
                TMatrix& A) override
        {
            this->refine(std::forward(A));
        }

        virtual TMatrix& refine(TMatrix& A) = 0;
};


template <class Matrix = DSparseMatrix<double>>
class CyclicPartitioner : public Partitioner<Matrix>
{
    public:

        CyclicPartitioner()
        {
            this->_procs = 1;
        };

        CyclicPartitioner(int procs)
        {
            this->_procs = procs;
        };

        virtual ~CyclicPartitioner() override {

        };

        virtual Matrix& partition(
                Matrix& A) override
        {
            using Image = typename Matrix::image_type;

            // repartition A
            A.setDistributionScheme(Partitioning::custom, this->_procs);

            vector<unique_ptr<Image>> new_images;
            for (int i = 0; i < this->_procs; ++i)
                new_images.push_back(std::make_unique<Image>());

            // FIXME, make this inplace for the first proc
            auto& images = A.getMutableImages();
            int idx = 0;
            for (auto& image : images)
                for(auto& trip : *image) {
                    new_images[idx++ % this->_procs]->pushTriplet(trip);
                }

            images.resize(this->_procs);
            for(int i = 0; i < this->_procs; ++i)
                images[i].reset(new_images[i].release());

            return A;
        }
};

}
