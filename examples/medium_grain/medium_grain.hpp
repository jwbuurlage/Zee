#include <zee.hpp>

#include <random>

#include <memory>
using std::shared_ptr;

template <class TMatrix = Zee::DSparseMatrix<double>>
class MGPartitioner : Zee::IterativePartitioner<TMatrix>
{
    public:
        virtual TMatrix& refine(TMatrix& A) override
        {
            return A;
        }
};
