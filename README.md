# Zee; A distributed matrix library and partitioning framework

<p align="center">
<img width="400px" alt="Zee is a matrix library and partitioning framework" src="https://raw.githubusercontent.com/jwbuurlage/zee/develop/docs/sphinx/img/zee.png" />
</p>

Sparse matrix partitioning and distributed numerical linear algebra.

Zee is a framework for implementing parallel algorithms for large-scale
linear algebra operations. The main focus is on facilitating general
partitioning methods for sparse matrices.

## Example usage:

```CPP
#include <zee.hpp>

// ...

DSparseMatrix<> A("matrix.mtx", 4);
DVector<> v{A.getCols(), 1.0};
DVector<> u{A.getRows(), 1.0};

u = A * v;

```

## License

Zee is released under the LGPLv3. See the file `COPYING` for details.

## Using Zee

Zee is a header-only `C++` library, and can be used by including `<zee.hpp>`
in your program. For detailed instructions we refer to the [documentation][dummy].

The `master` branch contains the latest release. An (unstable) snapshot of the current development can be found in the `develop` branch.

## Dependencies and requirements

The library has been tested on Linux using recent versions of Clang and GCC with support for C++14. For multithreading support the *pthreads* library is required.

Below is a list of optional dependencies:

- *Python* and the Python library *matplotlib* are used for scripts, in particular for plotting.
- [*Catch*][catch] is used for the unit tests.

[catch]: https://github.com/philsquared/Catch
[dummy]: #

## Contributing

See the file `contributing/CONTRIBUTING.md` for details on how to contribute.
