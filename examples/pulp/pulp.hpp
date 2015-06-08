#include <zee.hpp>

template <class TMatrix = Zee::DSparseMatrix<double>>
class PulpPartitioner : Zee::Partitioner<TMatrix>
{
    public:
        virtual TMatrix& refine(TMatrix& A)
        {
            auto& images = A.getMutableImages();
            for(auto& image : images) {
                // FIXME we need some kind of 'compute unit' associated with
                // a unit, for now we simulate this here (using threads?).

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
                //
                // TODO: Maybe first work on graph structure, write
                // something in tex about graph vs hypergraph partitioning
                // etc.

                // 1. Choose a random vertex to give away
                // 2. Ask all processors to give us their count
                // 3. Give it away to the major processor
                //
                // can we cache this information while spmving?

                // FIXME: rough sketch of code
                auto trip = image.getElement(rand() % image.size());
                auto col = trip.getCol();
                auto row = trip.getRow();

                // put in some kind of request
                vector<int> counts(A.getProcs() - 1, 0);
                auto countLambda = [] (TMatrixImage _img) {
                    auto count = 0; // FIXME: TIdx or auto
                    for (auto& a : _img) {
                        if (a.getCol() == col || a.getRow() == row)
                            count++;
                    }
                    return count;
                };

                A.compute(countLambda, counts);

                auto max = counts.maxIndex();
                images.yield(trip, image);
            }
        }
}
