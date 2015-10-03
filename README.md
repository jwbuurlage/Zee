Zee
===

Sparse matrix partitioning and distributed numerical linear algebra.

Zee is a framework for implementing parallel algorithms for large-scale
linear algebra operations. The main focus is on facilitating general
partitioning methods for sparse matrices.

Building
========

Zee is a header only `C++` library, and can be used by including `<zee.hpp>`
in your program. For detailed instructions we refer to the [documentation][dummy].
For multithreading support the `pthreads` library is required on linux.

Dependencies and requirements
=============================

Zee has only been tested on Linux. Below is a list of (optional) dependencies:

- `python` and the python package `matplotlib` are used for plotting.
- [`Catch`][dummy] is used for unit tests and is bundled with Zee. (legal?)

[dummy]: #
