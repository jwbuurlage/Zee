# General improvements
[ ] Write virtual container with iterators for triplets generation
[ ] Make triplet iteration more general, is it necessary to add iteration
    in template argument?

# Parallelilization
[ ] Formalize worker threads through some common provider interface.
    
# SpMV
[ ] Add multiple parallelization providers and SpMD implementations, make sure
    it works on general systems (not only multi-core systems, or via pthreads)
    In particular implement for BSP.

# Distributed types
[ ] Generalize base class even further, add support for distribution at highest
    possible level.
[ ] Generalize (sparse) iteration, make inheritance scheme.

# Partitioning
[ ] Partioners should modify or clone existing distributed structures.
    Need to devise an interface that makes partitioning as general as possible.
[ ] Consider Iterative Refinement (IR) methods and how to implement them.

# Matrices & Testing
[ ] Support for e.g. Matrix Market format
[ ] Benchmarking and testing, model communication volume and load imbalance.
