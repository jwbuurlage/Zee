#include <zee.hpp>
#include "partitioners/pulp.hpp"

#include "args/argparse.hpp"
#include "report/report.hpp"

using namespace Zee;

int main(int argc, char* argv[])
{
    using TVal = float;
    using TIdx = uint64_t;

    //----------------------------------------------------------------------------
    // parse cli args

    auto args = ArgParse(argc, argv);

    auto matrices = args.asList("--matrices");
    ZeeLogVar(matrices);

    HGModel model = HGModel::fine_grain;
    auto modelName = args.asString("--hg");
    if (modelName == "fine_grain")
        model = HGModel::fine_grain;
    else if (modelName == "row_net")
        model = HGModel::row_net;
    else if (modelName == "column_net")
        model = HGModel::column_net;
    else {
        ZeeLogInfo << "Using default hypergraph model (fine-grain)" << endLog;
    }

    TIdx procs = args.as<TIdx>("--procs");
    if (procs == 0) {
        procs = 2;
        ZeeLogInfo << "Using default number of processors (" << procs << ")"
                   << endLog;
    }

    auto epsilon = args.as<double>("--eps");
    if (epsilon == 0) {
        epsilon = 0.03;
        ZeeLogInfo << "Using default tolerance: " << epsilon << endLog;
    }

    bool plot = args.wasPassed("--plot");
    bool benchmark = args.wasPassed("--benchmark");

    TIdx runs = args.as<TIdx>("--runs");
    if (runs == 0) {
        runs = 1;
        ZeeLogInfo << "Using default number of runs (" << runs << ")" << endLog;
    }

    //----------------------------------------------------------------------------

    auto report = Report("Hyper-PuLP", "matrix");
    report.addColumn("m");
    report.addColumn("n");
    report.addColumn("N");
    report.addColumn("V_HP");
    report.addColumn("V_C");

    auto bench = Zee::Benchmark("PuLP");
    for (auto matrix : matrices) {
        report.addRow(matrix);
        bench.phase(matrix);

        // Initialize the matrix from file
        DSparseMatrix<TVal, TIdx> A{
            "data/matrices/" + matrix  + ".mtx",
            1
        };

        report.addResult(matrix, "m", A.getRows());
        report.addResult(matrix, "n", A.getCols());
        report.addResult(matrix, "N", A.nonZeros());

        //ZeeAssert(A.getRows() == A.getCols());

        auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
        auto u = DVector<TVal, TIdx>{A.getRows()};

        PulpPartitioner<decltype(A)> pA(A, procs, epsilon);
        pA.initialize(A, model);
        pA.initialPartitioning(5);

        std::vector<TIdx> communicationVolumes;
            communicationVolumes.push_back(A.communicationVolume());
        for (int i = 1; i <= 1000; ++i) {
            pA.refineWithIterations(A.nonZeros());
            auto comVol = A.communicationVolume();
            communicationVolumes.push_back(comVol);
            if (communicationVolumes[i] == communicationVolumes[i - 1])
                break;
        }
        ZeeLogVar(communicationVolumes);
        report.addResult(matrix, "V_HP", communicationVolumes.back());

        //A.spy("steam", true);
        if (plot) {
            auto p = Zee::Plotter<TIdx>();
            p["xlabel"] = "iterations";
            p["ylabel"] = "$V$";
            p["title"] = "Communication Volume PuLP ";
            p.addLine(communicationVolumes, "communicationVolumes");
            p.plot("comvol_pulp_" + matrix, true);
        }

        ZeeLogVar(A.loadImbalance());

        CyclicPartitioner<decltype(A)> cyclic(procs);
        cyclic.partition(A);

        report.addResult(matrix, "V_C", A.communicationVolume());

    }
    if (!benchmark) {
        bench.silence();
    }
    bench.finish();

    report.print();

//    for (auto matrix : matrices) {
//        // Initialize the matrix from file
//        DSparseMatrix<TVal, TIdx> A{
//            "data/matrices/" + matrix  + ".mtx",
//            1
//        };
//
//        MGPartitioner<decltype(A)> mgPartitioner;
//        mgPartitioner.partition(A);
//        A.spy(matrix + "_mg");
//
//        ZeeLogInfo << "MG: \t" << A.communicationVolume() << endLog;
//        while (!mgPartitioner.locallyOptimal()) {
//            mgPartitioner.refine(A);
//        }
//        ZeeLogVar(A.communicationVolume());
//    }

//    GreedyVectorPartitioner<decltype(A), decltype(v)> pVecs(A, v, u);
//    pVecs.partition();
//    pVecs.localizeMatrix();
//
//    u = A * v;

    return 0;
}
