#include <zee.hpp>
#include "partitioners/pulp.hpp"

using namespace Zee;

int main()
{
    using TVal = float;
    using TIdx = uint32_t;
    std::string matrix = "fpga_dcop_05";
    TIdx procs = 2;

    double epsilon = 0.03;

    // Initialize the matrix from file
    DSparseMatrix<TVal, TIdx> A{
        "data/matrices/" + matrix  + ".mtx",
        1
    };

    ZeeLogInfo << "Matrix of size (" << A.getRows() << " x " << A.getCols()
               << ")" << endLog;

    //ZeeAssert(A.getRows() == A.getCols());

    auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
    auto u = DVector<TVal, TIdx>{A.getRows()};

    PulpPartitioner<decltype(A)> pA(A, procs, epsilon);
    pA.initialize(A, HGModel::fine_grain);

    std::vector<TIdx> communicationVolumes;
        communicationVolumes.push_back(A.communicationVolume());
    for (int i = 1; i <= 1000; ++i) {
        pA.refineWithIterations(A.nonZeros());
        auto comVol = A.communicationVolume();
        ZeeLogVar(comVol);
        communicationVolumes.push_back(comVol);
        if (communicationVolumes[i] == communicationVolumes[i - 1])
            break;
    }

    //A.spy("steam", true);
    auto p = Zee::Plotter<TIdx>();
    p["xlabel"] = "iterations";
    p["ylabel"] = "$V$";
    p["title"] = "Communication Volume PuLP ";
    p.addLine(communicationVolumes, "communicationVolumes");
    p.plot("comvol_pulp_" + matrix);

    ZeeLogVar(communicationVolumes);
    ZeeLogVar(A.loadImbalance());

    ZeeLogInfo << "Cyclic for comparison:" << endLog;
    CyclicPartitioner<decltype(A)> cyclic(procs);
    cyclic.partition(A);
    ZeeLogVar(A.communicationVolume());
    ZeeLogVar(A.loadImbalance());

//    GreedyVectorPartitioner<decltype(A), decltype(v)> pVecs(A, v, u);
//    pVecs.partition();
//    pVecs.localizeMatrix();
//
//    u = A * v;

    return 0;
}
