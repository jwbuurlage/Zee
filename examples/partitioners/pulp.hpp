#include <zee.hpp>
#include <random>

// FIXME: TMP
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <memory>
using std::shared_ptr;

template <class TMatrix = Zee::DSparseMatrix<double>>
class PulpPartitioner : Zee::IterativePartitioner<TMatrix>
{
    public:
        virtual TMatrix& refine(TMatrix& A) override
        {
            std::random_device rd;
            std::mt19937 mt(rd());
            std::uniform_int_distribution<typename TMatrix::index_type>
                randmil(0, 1'000'000);

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

                auto size = image->nonZeros();
                if (size == 0) {
                    cerr << "ERROR: Matrix image empty!" << endl;
                    continue;
                }
                auto element_index = randmil(mt) % size; 

                // FIXME: rough sketch of code
                auto trip = image->getElement(element_index);
                auto col = trip.col();
                auto row = trip.row();

                // put in some kind of request
                auto countLambda = [col_ = col, row_ = row]
                    (std::shared_ptr<typename TMatrix::image_type> _img)
                    -> typename TMatrix::index_type
                {
                    typename TMatrix::index_type count = 0;
                    for (auto& a : *_img) {
                        if (a.col() == col_ || a.row() == row_)
                            count++;
                    }
                    return count;
                };

                auto counts = A.template compute<
                    typename TMatrix::index_type>(countLambda);

                auto max = 0;
                auto max_index = -1;
                for (typename TMatrix::index_type i = 0;
                        i < counts.size(); ++i) {
                    if (counts[i] > max) {
                        max = counts[i];
                        max_index = i;
                    }
                }

                // give element to max_index
                image->popElement(element_index);
                images[max_index]->pushTriplet(trip);
            }

            return A;
        }
};