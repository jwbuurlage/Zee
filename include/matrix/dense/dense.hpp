/*
File: include/matrix/dense/dense.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include <cstdint>

#include <vector>
#include <random>
#include <memory>
#include <atomic>
#include <mutex>
#include <map>

#include "../base/base.hpp"
#include "../sparse/sparse.hpp"
#include "../../common.hpp"
#include "../../operations/operations.hpp"
#include "../../default_types.hpp"

namespace Zee {

//-----------------------------------------------------------------------------
// VECTOR

template <typename Derived,
         typename TVal,
         typename TIdx>
class DVectorBase
    : public DMatrixBase<Derived, TVal, TIdx>
{
    public:
        using Base = DMatrixBase<Derived, TVal, TIdx>;
        using Base::operator=;

        DVectorBase(TIdx n, TVal defaultValue = 0)
            : Base(n, 1)
        { }
};

// FIXME: should be a specialization of a general dense matrix
// FIXME: saved as pairs? or just owners distributed cyclically
template <typename TVal = default_scalar_type,
         typename TIdx = default_index_type>
class DVector
    : public DVectorBase<DVector<TVal, TIdx>,
        TVal, TIdx>
{
    using Base = DVectorBase<DVector<TVal, TIdx>, TVal, TIdx>;

    public:
        using value_type = TVal;
        using index_type = TIdx;

        DVector(TIdx n, TVal defaultValue = 0)
            : Base(n, 1)
        {
            elements_.resize(n);
            std::fill(elements_.begin(), elements_.end(), defaultValue);
        }

        DVector(const DVector& other) :
            Base(other.size(), 1)
        {
            elements_ = other.elements_;
        }

        DVector(DVector&& other) :
            Base(other.size(), 1)
        {
            elements_ = std::move(other.elements_);
        }

        using Base::operator=;

        void operator= (DVector&& other) {
            this->rows_ = other.size();
            elements_ = std::move(other.elements_);
        }

        /** Return the size of the matrix (i.e. the number of rows) */
        TIdx size() const {
            return this->rows_;
        }

        TVal& operator[] (TIdx i)
        {
            return elements_[i];
        }

        const TVal& operator[] (TIdx i) const
        {
            return elements_[i];
        }

        /** The dot product */
        // FIXME parallellize
        TVal dot(const DVector<TVal, TIdx>& rhs) const {
            ZeeAssert(rhs.size() == size());

            const auto& lhs = *this;

            TVal sum = 0;
            for (TIdx i = 0; i < this->size(); ++i) {
                sum += lhs[i] * rhs[i];
            }

            return sum;
        }

        /** The (Euclidian) vector norm */
        TVal norm() const {
            TVal sum = 0;
            for (auto& element : elements_) {
                sum += element * element;
            }
            return sqrt(sum);
        }

        void reset() {
            for (auto& elem : elements_) {
                elem = 0;
            }
        }

        // Operator overloads and algorithm implementations

    private:
        std::vector<TVal> elements_;
};

// We add an operator such that we can log vectors
template <typename TVal, typename TIdx>
Logger& operator <<(Logger& lhs, const DVector<TVal, TIdx>& rhs) {
    auto sep = "";
    lhs << "[";
    for (TIdx i = 0; i < rhs.size(); ++i) {
        lhs << std::fixed << std::setprecision(2) << sep << rhs[i];
        sep = ", ";
    }
    lhs << "]";

    return lhs;
}

//-----------------------------------------------------------------------------
// MATRIX

// The dense base matrix
template <typename Derived,
         typename TVal = default_scalar_type,
         typename TIdx = default_index_type>
class DDenseMatrixBase
    : public DMatrixBase<Derived, TVal, TIdx>
{
    public:
        using Base = DMatrixBase<Derived, TVal, TIdx>;
        using Base::operator=;

        DDenseMatrixBase() { }

        DDenseMatrixBase(TIdx rows, TIdx cols)
            : Base(rows, cols)
        { }

        virtual TVal& at(TIdx i, TIdx j) = 0;
        virtual const TVal& at(TIdx i, TIdx j) const = 0;

        /** Obtain a spy image of the dense matrix */
        void spy(std::string title = "anonymous", bool show = false)
        {
            using std::endl;

            std::stringstream ss;
            ss << "data/spies/" << title << ".mtx";
            auto filename = ss.str();
            int i = 1;
            while(fileExists(filename)) {
                ss.str("");
                ss.clear();
                ss << "data/spies/" << title << "_" << i++ << ".mtx";
                filename = ss.str();
            }
            std::ofstream fout(filename);

            fout << "%%MatrixMarket matrix array real general" << endl;

            fout << this->getRows() << " " << this->getCols() << endl;

            for (TIdx j = 0; j < this->getCols(); ++j) {
                for (TIdx i = 0; i < this->getRows(); ++i) {
                    fout << this->at(i, j) << endl;
                }
            }

            ZeeLogInfo << "Spy saved to file: " << filename << Logger::end();

            if (show) {
                auto command = "./script/plot.py --showfile " + filename;
                std::system(command.c_str());
            }
        }

};

/* Contrary to a sparse matrix, we choose to not have a fixed distribution for
 * dense matrices, instead the distribution will be chosen in the implementation
 * of the algorithms.
 * When running a decentralized algorithm this class must be specialized. */
template <typename TVal = default_scalar_type,
         typename TIdx = default_index_type>
class DMatrix :
    public DDenseMatrixBase<DMatrix<TVal, TIdx>, TVal, TIdx>
{
    public:
        using Base = DDenseMatrixBase<DMatrix<TVal, TIdx>, TVal, TIdx>;
        using Base::operator=;

        DMatrix(TIdx rows, TIdx cols)
            : Base(rows, cols) {
            resize(rows, cols);
        }

        explicit DMatrix(std::string file)
            : Base(0, 0)
        {
            matrix_market::load(file, *this);
        }

        DMatrix(DMatrix& other)
        {
            *this = other;
        }

        DMatrix(DMatrix&& other)
        {
            *this = other;
        }

        void operator= (DMatrix& other) {
            elements_ = other.elements_;
            this->cols_ = other.getCols();
            this->rows_ = other.getRows();
        }

        void operator= (DMatrix&& other) {
            elements_ = std::move(other.elements_);
        }

        TVal& at(TIdx i, TIdx j)
        {
            ZeeAssert(i < this->rows_);
            ZeeAssert(j < this->cols_);
            return elements_[i][j];
        }

        const TVal& at(TIdx i, TIdx j) const
        {
            ZeeAssert(i < this->rows_);
            ZeeAssert(j < this->cols_);
            return elements_[i][j];
        }

        // column major order
        template<typename TInputIterator>
        void setFromValues(
            const TInputIterator& begin,
            const TInputIterator& end)
        {
            TInputIterator it = begin;
            for (TIdx j = 0; j < this->cols_; ++j) {
                for (TIdx i = 0; i < this->rows_; ++i)
                    elements_[i][j] = *(it++);
                }

            ZeeAssert(it == end);
        }

        void resize(TIdx rows, TIdx cols) override
        {

            Base::resize(rows, cols);

            elements_.resize(rows);
            for (auto& row : elements_) {
                row.resize(cols);
            }
        }

        void transpose()
        {
            // TODO: think about whether to explicitely (as it is now)
            // or implicitely transpose by simply using a bool flag and
            // changing dynamically at the various interfaces
            std::vector<std::vector<TVal>> new_elements;

            new_elements.resize(this->cols_);
            for (auto& col : new_elements) {
                col.resize(this->rows_);
            }

            Base::resize(this->cols_, this->rows_);

            for (TIdx i = 0; i < this->rows_; ++i) {
                for (TIdx j = 0; j < this->cols_; ++j) {
                    new_elements[i][j] = elements_[j][i];
                }
            }

            elements_ = std::move(new_elements);
        }

    private:
        // stored column major
        std::vector<std::vector<TVal>> elements_;
};

// We add an operator such that we can log dense matrices
template <typename TVal, typename TIdx>
Logger& operator<<(Logger& lhs, const DMatrix<TVal, TIdx>& rhs) {
    lhs << "\n";
    for (TIdx i = 0; i < rhs.getRows(); ++i) {
        lhs << "|";
        auto sep = "";
        for (TIdx j = 0; j < rhs.getCols(); ++j) {
            lhs << std::fixed << std::setprecision(2) << sep << rhs.at(i, j);
            sep = ", ";
        }
        lhs << "|";
        if (i != rhs.getRows() - 1)
            lhs << "\n";
    }
    return lhs;
}

// Operations involving dense vectors and matrices
#include "dense_operations.hpp"
#include "dense_matrix_operations.hpp"

} // namespace Zee
