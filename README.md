# Zee

### **A distributed matrix library and partitioning framework**

<p align="center">
<img width="400px" alt="Zee is a matrix library and partitioning framework" src="https://raw.githubusercontent.com/jwbuurlage/zee/develop/docs/sphinx/img/zee.png" />
</p>

**Zee** is a new partitioning framework. It is written in modern C++, and can be viewed as a unified software library for applications (initially linear algebra and in particular linear solvers) that also automatically and autonomously optimizes the running time of the operations involved on parallel systems, by means of balancing and partitioning.

The aim of Zee is to accomplish the following features:

* **Fully distributed**. Many partitioners are sequential, and are therefore not only unable to handle matrices that must be distributed from the start because of their size, they do not scale well either. A major goal of the framework is to be able to eventually do all computations in a completely distributed way, without needing centralized or shared memory.

* **Modular and extensible**. The library is built generally enough to support new partitioning methods as they are invented, without requiring to make major changes to the software. Furthermore, it is easy to extend its capabilities by adding additional partitioning methods or functionality, also for third-party users.

* **Maintainable and future-proof**. C++ has been a very successful language for performance-critical applications. The latest revisions of C++, C++11 and C++14 add many features to the language that make it particularly well-suited for the development of this framework.

* **Cross-platform**. It is designed to work on embedded systems, personal computers, and we plan to support distributed computing using e.g. MPI in a future release. In particular the library has been successfully tested on multi-core computers as well as the Parallella.

* **Integration** in applications. We focus on methods for iterative refinement of partitionings. The Zee framework is meant to be integrated in software employing these techniques. While first we target sparse matrix applications, in the future we could support more general applications using the same techniques.

The syntax and overall structure of Zee is based on [Eigen](http://eigen.tuxfamily.org), which is a 'C++ template library for linear algebra', so that it will feel familiar to use our new framework for developers familiar with this popular linear algebra library.

## Example usage:

Initializing vectors and matrices:

```CPP
#include <zee.hpp>

// ...

// matrices are initialized with a single image containing
// all zeros by default
auto A = DSparseMatrix<unsigned int, float>(rows, cols);

// we can construct sparse matrices by supplying triplets
std::vector<Triplet<unsigned int, float>> triplets;
A.setFromTriplets(triplets.begin(), triplets.end());

// they can also be loaded from a matrix market file
// by default the matrices are distributed cyclically over the
// `procs` processors
auto A = DSparseMatrix<unsigned int, float>("matrix.mtx", procs);

// dense matrices can be loaded in similarly
auto B = DMatrix<unsigned int, float>(rows, cols);
B.at(0, 0) = 1.0f;
auto B = DMatrix<unsigned int, float>("dense_matrix.mtx", procs);

// vectors are zero initialized by default, and can be
// viewed as (1 x n) dense matrices
auto v = DVector<unsigned int, float>(size);

// they can also be filled with a supplied value
auto value = 1.0f;
auto v = DVector<unsigned int, float>(size, value);
```

Linear algebra operations:

```CPP
// we initialize some objects
DVector<> v, w;
DSparseMatrix<> A;
DMatrix<> B;
DMatrix<> C;

// we can add and subtract vectors
DVector<> u1 = v + w;
DVector<> u2 = v - w;

// we can multiply matrices with vectors
DVector<> u1 = A * v;

// we can multiply matrices with matrices
DMatrix<> D = B * C;

// we can combine any number of operations as we see fit
DVector<> u1 = 2.0f * A * (v + w) + v;

// or even reuse existing vectors
w = A * v;
```

Partitioning matrices:
```CPP
// we initialize some objects
DVector<> v, w;
DSparseMatrix<> A, B;

// we can (bi-)partition the matrix A using e.g. medium-grain
MGPartitioner<decltype(A)> mediumGrain(epsilon);
mediumGrain.partition(A);

// we can refine the partitioner with MG-IR
while (!mediumGrain.locallyOptimal()) {
    mediumGrain.refine(A);
}

// we can also use Hyper-PuLP to partition
PulpPartitioner<decltype(B)> pulp(B, procs, epsilon);
pulp.initialize(HGModel::row_net);

// we obtain an initial partitioning
pulp.initialPartitioning(iterations);

// and now we can perform refinement iterations using
pulp.refine();

// if we want to perform a SpMV we need to partition the vector
GreedyVectorPartitioner<decltype(A), decltype(b)> greedy(A, v, w);
greedy.partition();

// we localize the indices of the matrix for efficiency
greedy.localizeMatrix();

// now we can perform the SpMV
v = A * w;
```

For more examples, see the [documentation](#).

## License

Zee is released under the LGPLv3. See the file `COPYING` for details.

## Using Zee

Zee is a header-only C++ library, and can be used by including `<zee.hpp>`
in your program. For detailed instructions we refer to the [documentation](#).

The `master` branch contains the latest release. An (unstable) snapshot of the current development can be found in the `develop` branch.

## Dependencies and requirements

The library has been tested on Linux using recent versions of Clang and GCC with support for C++14. For multithreading support the *pthreads* library is required.

Below is a list of optional dependencies:

- *Python* and the Python library *matplotlib* are used for scripts, in particular for plotting.
- [*Catch*](https://github.com/philsquared/Catch) is used for the unit tests.

## Contributing

See the file `contributing/CONTRIBUTING.md` for details on how to contribute.
