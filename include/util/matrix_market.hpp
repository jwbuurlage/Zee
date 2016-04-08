/*
File: include/matrix_market.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

// # TODO
// [ ] Support the extended format

#pragma once

#include <fstream>
#include <string>
#include <sstream>

#include "jw.hpp"

#include "matrix/base/base.hpp"

namespace Zee {

namespace matrix_market
{
    enum info : int {
        coordinate     = 1 << 0,
        array          = 1 << 1,
        real           = 1 << 2,
        pattern        = 1 << 3,
        general        = 1 << 4,
        symmetric      = 1 << 5,
        skew_symmetric = 1 << 6
    };

    template <typename Derived, typename TVal, typename TIdx>
    void loadMatrix(int info, std::ifstream& fs, DSparseMatrixBase<Derived, TVal, TIdx>& target)
    {
        if (info & matrix_market::info::array) {
            JWLogError << "Trying to load dense matrix .mtx "
                "file format into a sparse matrix." << endLog;
        }

        std::string line;

        while (!fs.eof()) {
            std::getline(fs, line);
            if (line[0] != '%')
                break;
        }

        TIdx M = 0;
        TIdx N = 0;
        TIdx L = 0;

        std::stringstream line_stream(line);

        // not a comment, if we have yet read M, N, L
        line_stream >> M >> N >> L;

        std::vector<Triplet<TVal, TIdx>> coefficients;

        // FIXME: why does this corrupt data
        //coefficients.reserve(L);

        // read matrix coordinates
        for (TIdx i = 0; i < L; ++i) {
            TIdx row, col;
            TVal value = (TVal)1;

            std::getline(fs, line);
            std::stringstream line_stream(line);
            line_stream >> row >> col;

            if (!(info & matrix_market::info::pattern)) {
                line_stream >> value;
            }

            coefficients.push_back(Triplet<TVal, TIdx>(row - 1, col - 1, value));
        }

        if (info & matrix_market::info::symmetric) {
            for (auto& trip : coefficients) {
                // we do not want to duplicate the diagonal
                if (trip.col() == trip.row()) {
                    continue;
                }

                coefficients.push_back(Triplet<TVal, TIdx>(
                    trip.col(), trip.row(), trip.value()));
            }
        } else if (info & matrix_market::info::skew_symmetric) {
            for (auto& trip : coefficients) {
                // we do not want to duplicate the diagonal
                if (trip.col() == trip.row())
                    continue;

                coefficients.push_back(Triplet<TVal, TIdx>(
                    trip.col(), trip.row(), -trip.value()));
            }
        }

        target.resize(M, N);
        target.setFromTriplets(coefficients.begin(), coefficients.end());
    }

    template <typename Derived, typename TVal, typename TIdx>
    void loadMatrix(int info, std::ifstream& fs, DMatrix<TVal, TIdx>& target)
    {

        if (info & matrix_market::info::coordinate) {
            JWLogError << "Trying to load sparse matrix .mtx "
                "file format into a dense matrix." << endLog;
        }

        std::string line;

        while (!fs.eof()) {
            std::getline(fs, line);
            if (line[0] != '%')
                break;
        }

        TIdx M = 0;
        TIdx N = 0;

        std::stringstream line_stream(line);

        // read M, N
        line_stream >> M >> N;

        std::vector<TVal> values;

        // read matrix coordinates
        for (TIdx i = 0; i < M * N; ++i) {
            TVal value = (TVal)1;

            std::getline(fs, line);
            std::stringstream line_stream(line);
            line_stream >> value;

            values.push_back(value);
        }

        target.resize(M, N);
        target.setFromValues(values.begin(), values.end());
    }

    /** Load matrix from .mtx format */
    template <typename Derived, typename TVal, typename TIdx>
    void load(std::string file, DMatrixBase<Derived, TVal, TIdx>& target)
    {
        JWLogInfo << "Loading matrix from file: " << file << endLog;

        std::ifstream fs(file);

        std::string header;
        std::getline(fs, header);
        std::stringstream header_stream(header);

        // first check if file is a legitimate MatrixMarket file
        std::string s, t;
        header_stream >> s >> t;

        if (s != "%%MatrixMarket" || t != "matrix") {
            JWLogError << "Not a valid MM file: " << file << endLog;
            exit(1);
        }

        // now follow a sequence of keywords
        int info = 0;
        while (!header_stream.eof()) {
            std::string keyword;
            header_stream >> keyword;

            if (s == "complex" ||
                    s == "Hermitian" ||
                    s == "integer") {
                JWLogError << "Unsupported keyword encountered: " << s << endLog;
                return;
            }

            if (keyword == "coordinate") {
                info |= matrix_market::info::coordinate;
            } else if (keyword == "array") {
                info |= matrix_market::info::array;
            } else if (keyword == "real") {
                info |= matrix_market::info::real;
            } else if (keyword == "pattern") {
                info |= matrix_market::info::pattern;
            } else if (keyword == "general") {
                info |= matrix_market::info::general;
            } else if (keyword == "symmetric") {
                info |= matrix_market::info::symmetric;
            } else if (keyword == "skew_symmetric") {
                info |= matrix_market::info::skew_symmetric;
            }
        }

        loadMatrix<Derived, TVal, TIdx>(info, fs, target.derived());
    }
} //namespace matrix_market

} // namespace Zee
