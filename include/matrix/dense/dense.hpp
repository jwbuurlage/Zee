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
#include <ostream>

#include <atomic>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <vector>

#include "../../operations/operations.hpp"
#include "../../util/common.hpp"
#include "../../util/default_types.hpp"
#include "../base/base.hpp"
#include "../sparse/sparse.hpp"

namespace Zee {

//-----------------------------------------------------------------------------
// VECTOR

template <typename Derived, typename TVal, typename TIdx>
class DVectorBase : public DMatrixBase<Derived, TVal, TIdx> {
   public:
    using Base = DMatrixBase<Derived, TVal, TIdx>;
    using Base::operator=;

    DVectorBase(TIdx n) : Base(n, 1) {}

    virtual void reassign(TIdx, TIdx) = 0;
};

// FIXME: should be a specialization of a general dense matrix
// FIXME: saved as pairs? or just owners distributed cyclically
template <typename TVal = default_scalar_type,
          typename TIdx = default_index_type>
class DVector : public DVectorBase<DVector<TVal, TIdx>, TVal, TIdx> {
    using Base = DVectorBase<DVector<TVal, TIdx>, TVal, TIdx>;

   public:
    using value_type = TVal;
    using index_type = TIdx;

    DVector(TIdx n, TVal defaultValue = 0) : Base(n) {
        elements_.resize(n);
        owners_.resize(n);
        std::fill(elements_.begin(), elements_.end(), defaultValue);
        std::fill(owners_.begin(), owners_.end(), (TIdx)0);
    }

    DVector(const DVector& other) : Base(other.size()) {
        elements_ = other.elements_;
    }

    DVector(DVector&& other) : Base(other.size()) {
        elements_ = std::move(other.elements_);
    }

    using Base::operator=;

    void operator=(DVector&& other) {
        this->rows_ = other.size();
        elements_ = std::move(other.elements_);
    }

    /** Return the size of the matrix (i.e. the number of rows) */
    TIdx size() const { return this->rows_; }

    TVal& operator[](TIdx i) { return elements_[i]; }

    const TVal& operator[](TIdx i) const { return elements_[i]; }

    /** The dot product */
    // FIXME parallellize
    TVal dot(const DVector<TVal, TIdx>& rhs) const {
        JWAssert(rhs.size() == size());

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

    void reassign(TIdx element, TIdx processorTarget) override {
        owners_[element] = processorTarget;
    }

    const std::vector<TIdx>& getOwners() const { return owners_; }

    bool operator==(const DVector<TVal, TIdx>& rhs) const {
        auto& lhs = (*this);
        if (lhs.size() != rhs.size()) return false;
        for (size_t i = 0; i < lhs.size(); ++i) {
            if (lhs[i] != rhs[i]) return false;
        }
        return true;
    }

   private:
    std::vector<TVal> elements_;
    std::vector<TIdx> owners_;
};

// We add an operator such that we can log vectors
template <typename TVal, typename TIdx>
jw::Logger& operator<<(jw::Logger& lhs, const DVector<TVal, TIdx>& rhs) {
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
template <typename Derived, typename TVal = default_scalar_type,
          typename TIdx = default_index_type>
class DDenseMatrixBase : public DMatrixBase<Derived, TVal, TIdx> {
   public:
    using Base = DMatrixBase<Derived, TVal, TIdx>;
    using Base::operator=;

    DDenseMatrixBase() {}

    DDenseMatrixBase(TIdx rows, TIdx cols) : Base(rows, cols) {}

    virtual TVal& at(TIdx i, TIdx j) = 0;
    virtual const TVal& at(TIdx i, TIdx j) const = 0;

    /** Obtain a spy image of the dense matrix */
    void spy(std::string title = "anonymous", bool show = false) {
        using std::endl;

        std::stringstream ss;
        ss << "data/spies/" << title << ".mtx";
        auto filename = ss.str();
        int i = 1;
        while (fileExists(filename)) {
            ss.str("");
            ss.clear();
            ss << "data/spies/" << title << "_" << i++ << ".mtx";
            filename = ss.str();
        }
        std::ofstream fout(filename);

        fout << "%%MatrixMarket matrix array real general" << endl;

        fout << this->getRows() << " " << this->getCols() << endl;

        for (TIdx j = 0; j < this->getCols(); ++j) {
            for (TIdx k = 0; k < this->getRows(); ++k) {
                fout << this->at(k, j) << endl;
            }
        }

        JWLogInfo << "Spy saved to file: " << filename << endLog;

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
class DMatrix : public DDenseMatrixBase<DMatrix<TVal, TIdx>, TVal, TIdx> {
   public:
    using Base = DDenseMatrixBase<DMatrix<TVal, TIdx>, TVal, TIdx>;
    using Base::operator=;

    DMatrix(TIdx rows, TIdx cols) : Base(rows, cols) { resize(rows, cols); }

    explicit DMatrix(std::string file) : Base(0, 0) {
        matrix_market::load(file, *this);
    }

    DMatrix(DMatrix& other) : Base(other.getRows(), other.getCols()) {
        *this = other;
    }

    DMatrix(DMatrix&& other) : Base(other.getRows(), other.getCols()) {
        *this = other;
    }

    void operator=(DMatrix& other) {
        elements_ = other.elements_;
        this->cols_ = other.getCols();
        this->rows_ = other.getRows();
    }

    void operator=(DMatrix&& other) { elements_ = std::move(other.elements_); }

    TVal& at(TIdx i, TIdx j) {
        JWAssert(i < this->rows_);
        JWAssert(j < this->cols_);
        return elements_[i][j];
    }

    const TVal& at(TIdx i, TIdx j) const {
        JWAssert(i < this->rows_);
        JWAssert(j < this->cols_);
        return elements_[i][j];
    }

    // row major order
    template <typename TInputIterator>
    void setFromValues(const TInputIterator& begin, const TInputIterator& end) {
        TInputIterator it = begin;
        for (TIdx i = 0; i < this->rows_; ++i)
            for (TIdx j = 0; j < this->cols_; ++j) elements_[i][j] = *(it++);

        JWAssert(it == end);
    }

    void resize(TIdx rows, TIdx cols) override {
        Base::resize(rows, cols);

        elements_.resize(rows);
        for (auto& row : elements_) {
            row.resize(cols);
        }
    }

    void transpose() {
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
std::ostream& operator<<(std::ostream& lhs, const DMatrix<TVal, TIdx>& rhs) {
    for (TIdx i = 0; i < rhs.getRows(); ++i) {
        lhs << "|";
        auto sep = "";
        for (TIdx j = 0; j < rhs.getCols(); ++j) {
            lhs << std::fixed << std::setprecision(2) << sep << rhs.at(i, j);
            sep = ", ";
        }
        lhs << "|";
        if (i != rhs.getRows() - 1) lhs << "\n";
    }
    return lhs;
}

// We add an operator such that we can log dense matrices
template <typename TVal, typename TIdx>
std::ostream& operator<<(std::ostream& lhs, const DVector<TVal, TIdx>& rhs) {
    lhs << "|";
    auto sep = "";
    for (TIdx j = 0; j < rhs.size(); ++j) {
        lhs << std::fixed << std::setprecision(2) << sep << rhs[j];
        sep = ", ";
    }
    lhs << "|";
    return lhs;
}

// Operations involving dense vectors and matrices
#include "dense_matrix_operations.hpp"
#include "dense_operations.hpp"

}  // namespace Zee
