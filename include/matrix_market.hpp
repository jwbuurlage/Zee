/*
File: include/sparse_matrix.h

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

#include "matrix/sparse.hpp"
#include "logging.hpp"

namespace Zee {

enum MMInfo : int {
    coordinate     = 1 << 0,
    real           = 1 << 1,
    pattern        = 1 << 2,
    general        = 1 << 3,
    symmetric      = 1 << 4,
    skew_symmetric = 1 << 5
};

/** Load a sparse matrix from MM */
template <typename TVal, typename TIdx>
DSparseMatrix<TVal, TIdx> fromMatrixMarket(std::string file, TIdx procs)
{
    auto n = 40;

    ZeeLogInfo << "Loading matrix from file: " << file << endLog;

    std::ifstream fs(file);

    std::string header;
    std::getline(fs, header);
    std::stringstream header_stream(header);

    // first check if file is a legitimate MatrixMarket file
    std::string s, t;
    header_stream >> s >> t;

    if (s != "%%MatrixMarket" || t != "matrix") {
        ZeeLogError << "Not a valid MM file: " << file << endLog;
        return DSparseMatrix<TVal, TIdx>(1, 1);
    }

    // now follow a sequence of keywords
    int info = 0;
    while (!header_stream.eof()) {
        std::string keyword;
        header_stream >> keyword;

        if (s == "array" ||
                s == "complex" ||
                s == "Hermitian" ||
                s == "integer") {
            ZeeLogError << "Unsupported keyword encountered: " << s << Logger::end();
            return DSparseMatrix<TVal, TIdx>(1, 1);
        }

        if (keyword == "coordinate") {
            info |= MMInfo::coordinate;
        } else if (keyword == "real") {
            info |= MMInfo::real;
        } else if (keyword == "pattern") {
            info |= MMInfo::pattern;
        } else if (keyword == "general") {
            info |= MMInfo::general;
        } else if (keyword == "symmetric") {
            info |= MMInfo::symmetric;
        } else if (keyword == "skew_symmetric") {
            info |= MMInfo::skew_symmetric;
        }

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

    // not a comment, if we have yet read N, M, L
    line_stream >> M >> N;
    if (line_stream.eof()) {
        // error dense matrix
        ZeeLogError << "Dense matrix format not supported" << Logger::end();
        return DSparseMatrix<TVal, TIdx>(1, 1);
    }

    line_stream >> L;

    vector<Triplet<TVal, TIdx>> coefficients;
    coefficients.reserve(n);

    // read matrix coordinates
    for (TIdx i = 0; i < L; ++i) {
        TIdx row, col;
        TVal value = (TVal)1;

        std::getline(fs, line);
        std::stringstream line_stream(line);
        line_stream >> row >> col;

        if (!(info & MMInfo::pattern)) {
            line_stream >> value;
        }

        coefficients.push_back(Triplet<TVal, TIdx>(row - 1, col - 1, value));
    }

    if (info & MMInfo::symmetric) {
        for (auto& trip : coefficients) {
            coefficients.push_back(Triplet<TVal, TIdx>(
                trip.col(), trip.row(), trip.value()));
        }
    } else if (info & MMInfo::skew_symmetric) {
        for (auto& trip : coefficients) {
            coefficients.push_back(Triplet<TVal, TIdx>(
                trip.col(), trip.row(), -trip.value()));
        }
    }

    DSparseMatrix<TVal, TIdx> A(M, N);
    A.setDistributionScheme(Partitioning::cyclic, procs);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

} // namespace Zee
