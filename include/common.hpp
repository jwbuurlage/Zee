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

#include <map>
#include <atomic>

namespace Zee {

using std::make_pair;
using std::atomic;

template <typename T>
class counted_set :
    public std::map<T, T>
{
    public:
        counted_set() : std::map<T, T>() { }

        void raise(T key) {
            if (this->find(key) != this->end()) {
                (*this)[key] += 1;
            } else {
                this->insert(make_pair(key, 1));
            }
        }

        void lower(T key) {
            if ((*this)[key] > 1) {
                (*this)[key] -= 1;
            } else {
                this->erase(key);
            }
        }
};

// When atomic is used in a nested vector, it is constructed in 2 phases if 
// I understand correctly. Since atomic integrals have their move/copy ctors
// deleted, this does not compile. This wrapper simply adds copy/move ctors, but
// these are *not* atomic, and thus not thread-safe. Use with caution.
template <typename T>
class atomic_wrapper
{
    public:
        atomic<T> _a;

        atomic_wrapper() : _a(0) { }


        atomic_wrapper(const std::atomic<T> &a)
            : _a(a.load()) {}

        atomic_wrapper(const atomic_wrapper &other)
            : _a(other._a.load()) {}

        atomic_wrapper &operator=(const atomic_wrapper &other)
        {
            _a.store(other._a.load());
        }
};

} // namespace Zee
