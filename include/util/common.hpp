/*
File: include/common.hpp

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
#include <memory>
#include <mutex>
#include <condition_variable>

#include <sys/stat.h>

namespace Zee {

using std::make_pair;
using std::unique_ptr;

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

        atomic_wrapper() : a(0) { }

        atomic_wrapper(const std::atomic<T>& a)
            : a(a.load()) {}

        atomic_wrapper(const atomic_wrapper& other)
            : a(other.a.load()) {}

        void operator=(const atomic_wrapper& other)
        {
            a.store(other.a.load());
        }

        void operator=(const T& other)
        {
            a = other;
        }

        std::atomic<T> a;
};

//* General factory class */
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

bool fileExists(std::string path)
{
    struct stat fileInfo;
    return stat(path.c_str(), &fileInfo) == 0;
}

// FIXME move to center
template <typename TIdx>
class Barrier {
    public:
        Barrier (TIdx procs = 0)
            : procs_(procs) {}

        inline void sync()
        {
            std::unique_lock<std::mutex> lock(mtx);
            count_++;
            if (count_ == procs_) {
                cv.notify_all();
                count_ = 0;
            } else {
                cv.wait(lock);
            }
        }

    private:
        std::mutex mtx;
        std::condition_variable cv;
        TIdx procs_ = 0;
        TIdx count_ = 0;
};

} // namespace Zee
