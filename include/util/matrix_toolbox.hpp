/*
This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/
#pragma once

#include "matrix/sparse/sparse.hpp"

namespace Zee {

//-----------------------------------------------------------------------------
// Convenience functions (MATLAB syntax)
//-----------------------------------------------------------------------------

/** Create identity matrix as sparse matrix */
template <typename TIdx>
DSparseMatrix<default_scalar_type, TIdx> eye(TIdx n, TIdx procs) {
    std::vector<Triplet<default_scalar_type, TIdx>> coefficients;
    coefficients.reserve(n);

    for (TIdx i = 0; i < n; ++i)
        coefficients.push_back(Triplet<default_scalar_type, TIdx>(i, i, 1));

    DSparseMatrix<default_scalar_type, TIdx> A(n, n);
    A.setDistributionScheme(partitioning_scheme::cyclic, procs);

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

/** Create a random sparse (n x m) matrix */
template <typename TIdx>
DSparseMatrix<double, TIdx> rand(TIdx m, TIdx n, TIdx procs, double density) {
    std::vector<Triplet<double, TIdx>> coefficients;
    coefficients.reserve((int)(n * m * density));

    double mu = 1.0 / density + 0.5;
    double sigma = 0.5 * mu;

    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::normal_distribution<double> gauss(mu, sigma);

    TIdx j = static_cast<int>(gauss(mt)) / 2;
    j = std::max(j, (TIdx)1);
    TIdx i = 0;
    while (true) {
        coefficients.push_back(
            Triplet<double, TIdx>(i, j, (1.0 + 10.0 * dist(mt))));

        int offset = static_cast<int>(gauss(mt));
        offset = std::max(offset, 1);
        j += offset;

        while (j >= n) {
            j = j - n;
            i++;
        }

        if (i >= m) break;
    }

    DSparseMatrix<double, TIdx> A(m, n);
    A.setDistributionScheme(partitioning_scheme::random, procs);

    A.setFromTriplets(coefficients.begin(), coefficients.end());

    return A;
}

}  // namespace Zee
