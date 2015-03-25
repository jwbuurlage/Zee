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

#include <iterator>
#include <type_traits>

#include <memory>
using std::unique_ptr;


namespace Zee
{

template <typename TVal, typename TIdx>
class Triplet;

/** Storage type for sparse matrix (image). */
enum StorageType
{
    S_TRIPLETS,
    S_COMPRESSED_ROW, // CRS
    S_COMPRESSED_COLUMN // CCS
};

template <typename TVal, typename TIdx>
class DSparseStorage
{
    public:
        DSparseStorage() { }
        virtual ~DSparseStorage() { }

        virtual void pushTriplet(Triplet<TVal, TIdx> t) { }

        /** We define iterators and constant iterators for getting out
          * triplets */
        typedef std::iterator<
                std::bidirectional_iterator_tag,
                Triplet<TVal, TIdx>> iterator;

        // We want to get triplets out of this
        virtual iterator& begin() = 0;
        virtual iterator& end() = 0;
};

template <typename TVal, typename TIdx>
class DSparseStorageTriplets : public DSparseStorage<TVal, TIdx>
{
    public:

        /** Base class for storage iteration */
        template <bool const_iter = true>
        class dual_iterator :
            public std::iterator<
                std::bidirectional_iterator_tag,
                Triplet<TVal, TIdx>> {
            public:
                // Define whether our data type is constant
                typedef typename std::conditional<const_iter,
                        const DSparseStorage<TVal, TIdx>*,
                        DSparseStorage<TVal, TIdx>*>::type StoragePointer;

                typedef typename std::conditional<const_iter,
                        const Triplet<TVal, TIdx>&,
                        Triplet<TVal, TIdx>&>::type TripletReference;

                /** Default constructor */
                dual_iterator(StoragePointer storage, TIdx i) :
                    _storage(storage),
                    _i(i) { }

                /** Copy constructor (const <-> regular conversion) */
                dual_iterator(const dual_iterator<false>& other) :
                    _storage(other._storage),
                    _i(other._i)  { }

                dual_iterator& operator--(int)
                {
                    // copy iterator, decrease this and return old copy
                    const dual_iterator old(*this);
                    --(*this);
                    return old;
                }

                dual_iterator& operator++(int)
                {
                    // copy iterator, increase this and return old copy
                    const dual_iterator old(*this);
                    ++(*this);
                    return old;
                }

                bool operator== (const dual_iterator& other) const
                {
                    return (_i == other._i);
                }

                /** Comparison operator for not equal
                  * @see operator==(const dual_iterator& other) const */
                bool operator!= (const dual_iterator& other) const
                {
                    return !(*this == other);
                }

                TripletReference operator*()
                {
                    return &_storage._triplets[_i];
                }

                dual_iterator& operator--() {
                    _i--;
                    return *this;
                }

                dual_iterator& operator++()
                {
                    _i++;
                    return *this;
                }

               /* now the copy constructor can access _storage for 
                * conversion */
               friend class dual_iterator<true>; 

            private:
                StoragePointer _storage;
                TIdx _i;
        };


        /** We define iterators and constant iterators for getting out
          * triplets */
        typedef dual_iterator<false> iterator_ts;
        typedef dual_iterator<true> const_iterator_ts;

        iterator_ts& begin()
        {
            // FIXME smart pointer???
            return *(new iterator_ts(this, 0));
        }

        iterator_ts& end()
        {
            // FIXME smart pointer???
            return *(new iterator_ts(this, _triplets.size()));
        }

        DSparseStorageTriplets() { };
        ~DSparseStorageTriplets() { };

        void pushTriplet(Triplet<TVal, TIdx> t)
        {
            // FIXME: does this get copied? want by ref, move?
            _triplets.push_back(t);
        }
        
    private:
        // FIXME: proper wasy to store this (reference?)
        vector<Triplet<TVal, TIdx>> _triplets;
};

}
