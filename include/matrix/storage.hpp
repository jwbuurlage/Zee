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

template <typename TVal, typename TIdx>
class DSparseStorageTriplets;

template <typename TVal, typename TIdx,
         class Derived = DSparseStorageTriplets<TVal, TIdx>>
class DSparseStorage;

/** FIXME: There should be a base class for storage iteration */
template <typename TVal, typename TIdx, bool const_iter = true>
class storage_iterator_triplets:
    public std::iterator<
        std::bidirectional_iterator_tag,
        Triplet<TVal, TIdx>> {
    public:
        // Define whether our data type is constant
        typedef typename std::conditional<const_iter,
                const DSparseStorageTriplets<TVal, TIdx>*,
                DSparseStorageTriplets<TVal, TIdx>*>::type
                    StoragePointer;

        typedef typename std::conditional<const_iter,
                const Triplet<TVal, TIdx>&,
                Triplet<TVal, TIdx>&>::type TripletReference;

        /** Default constructor */
        storage_iterator_triplets(StoragePointer storage, TIdx i) :
            _storage(storage),
            _i(i) { }

        /** Copy constructor (const <-> regular conversion) */
        storage_iterator_triplets(
                const storage_iterator_triplets<TVal, TIdx, false>& other) :
            _storage(other._storage),
            _i(other._i)  { }

        storage_iterator_triplets& operator--(int)
        {
            // copy iterator, decrease this and return old copy
            const storage_iterator_triplets old(*this);
            --(*this);
            return old;
        }

        storage_iterator_triplets& operator++(int)
        {
            // copy iterator, increase this and return old copy
            const storage_iterator_triplets old(*this);
            ++(*this);
            return old;
        }

        bool operator== (const storage_iterator_triplets& other) const
        {
            return (_i == other._i);
        }

        /** Comparison operator for not equal
          * @see operator==(const dual_iterator& other) const */
        bool operator!= (const storage_iterator_triplets& other) const
        {
            return !(*this == other);
        }

        TripletReference operator*()
        {
            return _storage->_triplets[_i];
        }

        storage_iterator_triplets& operator--() {
            _i--;
            return *this;
        }

        storage_iterator_triplets& operator++()
        {
            _i++;
            return *this;
        }

       /* now the copy constructor can access _storage for 
        * conversion */
       friend class storage_iterator_triplets<TVal, TIdx, true>; 

    private:
        StoragePointer _storage;
        TIdx _i;
};


template <typename TVal, typename TIdx, class TIterator>
class DSparseStorage
{
    public:
        DSparseStorage() { }
        virtual ~DSparseStorage() { }

        virtual void pushTriplet(Triplet<TVal, TIdx> t) { }

        /** We define iterators and constant iterators for getting out
          * triplets */
        virtual TIterator begin() = 0;
        virtual TIterator end() = 0;
};

template <typename TVal, typename TIdx>
class DSparseStorageTriplets :
    public DSparseStorage<TVal, TIdx,
        storage_iterator_triplets<TVal, TIdx, false>>
{
    public:
        /** We define iterators and constant iterators for getting out
          * triplets */
        typedef storage_iterator_triplets<TVal, TIdx, false> iterator;
        typedef storage_iterator_triplets<TVal, TIdx, true> const_iterator;

        iterator begin()
        {
            return iterator(this, 0);
        }

        iterator end()
        {
            return iterator(this, _triplets.size());
        }

        DSparseStorageTriplets() { }
        ~DSparseStorageTriplets() { }

        void pushTriplet(Triplet<TVal, TIdx> t)
        {
            // FIXME: does this get copied? want by ref, move?
            _triplets.push_back(t);
        }
        
        friend iterator;
        friend const_iterator;

    private:
        // FIXME: proper wasy to store this (reference?)
        vector<Triplet<TVal, TIdx>> _triplets;
};

}
