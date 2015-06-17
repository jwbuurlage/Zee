#include <zee.hpp>

using namespace Zee;

int main()
{
    DSparseMatrix<double, int> A = fromMM<double, int>("data/karate.mtx");

    A.spy();
    cout << A.loadImbalance() << endl;
    cout << A.communicationVolume() << endl;

    return 0;
}
