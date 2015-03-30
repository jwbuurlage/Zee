# Distributed types
[!] Generalize base class even further, add support for distribution at highest
possible level. In particular make (dense) vectors distributed and add
dense matrix type.
[!] Generalize (sparse) iteration, make inheritance scheme.

# Partitioning
[ ] Partioners should modify or clone existing distributed structures.
Need to devise an interface that makes partitioning as general as possible.
[ ] Consider Iterative Refinement (IR) methods and how to implement them.

# Parallelilization
[ ] Formalize worker threads through some common provider interface.
[ ] Formalize (parallel) operation generalization. How to specify for different
(parallel) providers
    
# SpMV
[ ] Add multiple parallelization providers and SpMD implementations, make sure
it works on general systems (not only multi-core systems, or via pthreads)
In particular implement for BSP.
[ ] IMPORTANT: completely parallel version of SPMD, without O(n) storage
requirements and low computational complexity, also ability to cache?

# General improvements
[ ] Write virtual container with iterators for triplets generation
[ ] Make triplet iteration more general, is it necessary to add iteration in
templatized argument? Rather have storage templatized, and alias or typedef
iteration in storage (necessary anyway).
[ ] Generalize spy for large matrices
[ ] Make operations encapsulated in objects, and only perform them at an
= operation on a distributed type. (expression templates)  Requires careful
design.
[x] move `using` aliases to Zee namespace
[ ] Make logging mechanism

# Matrices & Testing
[ ] Support for e.g. Matrix Market format
[ ] Benchmarking `<chrono>` and (unit) testing,
[ ] Model communication volume and load imbalance.
[-] Output performance graphs, spy-ish distribution visualization
