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

#include <fstream>
#include <string>
#include <sstream>

#include "matrix/sparse.hpp"
#include "logging.hpp"

namespace Zee {

/** Create identity matrix as sparse matrix */
template <typename TVal, typename TIdx>
DSparseMatrix<TVal, TIdx> fromMM(std::string file)
{
    auto n = 40;
    auto procs = 4;

    logInfo("Loading matrix");

    std::ifstream fs(file);

    std::string header;
    std::getline(fs, header);
    std::stringstream header_stream(header);

    // first check if file is a legitimate MatrixMarket file
    std::string s, t;
    header_stream >> s >> t;

    if (s != "%%MatrixMarket" || t != "matrix") {
        logError("Not a valid MM file");
        return DSparseMatrix<TVal, TIdx>(1, 1);
    }

    // now follow a sequence of keywords

    line_stream >> s >> t >> u;

    if (s != "coordinate" || t != "real" || u != "general") {
        logError("Currently only real matrices in general coordinate "
                "representation are supported");
        return DSparseMatrix<TVal, TIdx>(1, 1);
    }

    logInfo(s);

    vector<Triplet<TVal, TIdx>> coefficients;
    coefficients.reserve(n);

    for (TIdx i = 0; i < n; ++i)
        coefficients.push_back(Triplet<TVal, TIdx>(i, i, 1.0));

    DSparseMatrix<TVal, TIdx> A(n, n);
    A.setDistributionScheme(Partitioning::cyclic, procs);
    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

} // namespace Zee
