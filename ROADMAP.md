TODAY:
------

[!] Benchmarking `<chrono>` 
[!] Output performance graphs
    - need to define filetypes for output data (see e.g. spy)
    o spy-ish distribution visualization
    o plots for communication volume
[!] Think of data output, write python scripts for plotting
[ ] Implement Multi-Level KL with row-net.


HIGH PRIORITY:
--------------

# Distributed types
[!] Generalize base class even further, add support for distribution at highest
possible level. In particular make (dense) vectors distributed and add
dense matrix type.
[ ] Add support for compressed storage

# Partitioning
[ ] Need to consider in-place (re)partitioning, reusing images in particular

# Images
[ ] Need to think of an approach to make an image something more physical.
In particular for distributed systems we require that the image is something
which exists for longer time spans and can be reused. We want to do this in
a portable way.
- Idea: serialize an image somehow, and extract it at a 'physical location'.
 * not portable (because we want to do this in a distributed fashion), unless
 * the storage is never really packed, but instead some kind of generator
   approach is used.

REFINEMENT:
-----------

# Parallelilization
[ ] Formalize worker threads through some common provider interface.
[ ] Formalize (parallel) operation generalization. How to specify for different
(parallel) providers

# SpMV
[ ] Add multiple parallelization providers and SpMD implementations, make sure
it works on general systems (not only multi-core systems, or via pthreads)
In particular implement for BSP.
[!] Completely parallel version of SPMD, without O(n) storage
requirements and low computational complexity, also ability to cache?

# General improvements
[ ] Generalize spy for large matrices
[ ] Make operations encapsulated in objects, and only perform them at an
= operation on a distributed type. (expression templates).  Requires careful
design.
[ ] Write virtual container with iterators for triplets generation
[ ] Underscores are ugly..

# (Test) Matrices & Testing / Benchmarking
[!] Unit testing, custom solution w/ Python
    - Unit tests using `script/test.py [category]`
    - Actual tests in  `test/category.c`
    - See unit tests in `munificent/wren`
    - OR use something like Catch, gtest, etc.

# Done
[x] move `using` aliases to Zee namespace
[x] Make triplet iteration more general, is it necessary to add iteration in
templatized argument? Rather have storage templatized, and alias or typedef
iteration in storage (necessary anyway).
[x] Generalize (sparse) iteration, make inheritance scheme. (solved by traits)
[x] Storage type in DSparseMatrix
[x] Partioners should modify or clone existing distributed structures.
Need to devise an interface that makes partitioning as general as possible.
[x] Images depend on matrix to be partitioned, can we use type aliases?
[x] Make Logging mechanism `logging.hpp`
[x] Support for e.g. Matrix Market format
[x] Model communication volume and load imbalance.
    - Do we explicitely construct (hyper)graph?
    - Possible to cache and (incrementally) update?
[x] Consider Iterative Refinement (IR) methods and how to implement them.
[x] Implement MG
