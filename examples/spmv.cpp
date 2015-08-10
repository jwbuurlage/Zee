#include <string>

#include <unpain_cpp.hpp>
#include <zee.hpp>

using namespace Zee;

int main()
{
    auto center = std::make_shared<Unpain::Center<int>>(2);
    std::string matrix = "steam3";

    ZeeLogInfo << "-- Starting SpMV example" << endLog;

    auto A = DSparseMatrix<double, int>(center, "data/matrices/" + matrix + ".mtx");
    A.spy();

    auto v = DVector<double, int>{center, A.getCols(), 1.0};

    auto u = DVector<double, int>{center, A.getRows()};

    u = A * v;

    ZeeLogVar(u[3]);

    return 0;
}
