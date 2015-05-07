#include <zee.hpp>

template <class TMatrix = Zee::DSparseMatrix<double>>
class PulpPartitioner : Zee::Partitioner<TMatrix>
{
    public:
        virtual TMatrix& refine(TMatrix& A)
        {
            auto& images = A.getMutableImages();
            for(auto& image : images) {
                // now we need to let every 'image' calculate
                // or refine its partitioning.
                //
                // what is required for this:
                // - Local knowledge of own elements
                // - Local knowledge of halo
                // - Local information about global structure(?)
                //    * as hints? p raw values?
                //
                // in particular we need to let every 'image' which
                // represents most generally a real processor with
                // a subset of the the distributed data, refine the
                // distribution locally
                //
                // This means we need to extend image to something
                // more physical.
                //
                // Let's write this out on paper.
            }
        }
}
