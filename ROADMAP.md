TODAY:
------

Partitioners that take a matrix, and (re)partition them.
- How do we use images, in particular -- storage!
- First generalize the iterators to some base iterator, and get rid of templates
- This is actually hard to do for abstract base classes....


HIGH PRIORITY:
--------------

# Distributed types
[!] Generalize base class even further, add support for distribution at highest
possible level. In particular make (dense) vectors distributed and add
dense matrix type.
[!] Generalize (sparse) iteration, make inheritance scheme.
[ ] Storage type in DSparseMatrix
[ ] Add support for compressed storage

# Partitioning
[ ] Partioners should modify or clone existing distributed structures.
Need to devise an interface that makes partitioning as general as possible.
[ ] Need to consider in-place (re)partitioning, reusing images in particular
[ ] Images depend on matrix to be partitioned, can we use type aliases?
[ ] Consider Iterative Refinement (IR) methods and how to implement them.

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
[!] Make Logging mechanism `logging.hpp`
[ ] Generalize spy for large matrices
[ ] Make operations encapsulated in objects, and only perform them at an
= operation on a distributed type. (expression templates).  Requires careful
design.
[ ] Write virtual container with iterators for triplets generation

# (Test) Matrices & Testing / Benchmarking
[ ] Support for e.g. Matrix Market format
[ ] Benchmarking `<chrono>` 
[ ] Unit testing, custom solution w/ Python
    - Unit tests using `script/test.py [category]`
    - Actual tests in  `test/category.c`
    - See unit tests in `munificent/wren`
    - OR use something like Catch, gtest, etc.
[ ] Model communication volume and load imbalance.
    - Do we explicitely construct (hyper)graph?
[-] Output performance graphs, spy-ish distribution visualization

# Done
[x] move `using` aliases to Zee namespace
[x] Make triplet iteration more general, is it necessary to add iteration in
templatized argument? Rather have storage templatized, and alias or typedef
iteration in storage (necessary anyway).
