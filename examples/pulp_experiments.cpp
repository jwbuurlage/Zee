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
    // parse CLI args

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
    bool randomize = args.wasPassed("--randomize");
    bool compareMg = args.wasPassed("--compare-mg");

    TIdx runs = args.as<TIdx>("--runs");
    if (runs == 0) {
        runs = 1;
        ZeeLogInfo << "Using default number of runs (" << runs << ")" << endLog;
    }

    TIdx iters = args.as<TIdx>("--iters");
    if (iters == 0) {
        iters = 1;
        ZeeLogInfo << "Using default number of iters (" << iters << ")" << endLog;
    }

    //----------------------------------------------------------------------------

    auto report = Report("Hyper-PuLP", "matrix");
    report.addColumn("m");
    report.addColumn("n");
    report.addColumn("N");
    report.addColumn("V_HP");
    report.addColumn("SD");
    report.addColumn("V_C");
    if (procs == 2)
        report.addColumn("V_MG");

    auto bench = Zee::Benchmark("PuLP");

    std::vector<double> improvements;
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

        std::vector<double> Vs;
        for (TIdx run = 0; run < runs; ++run) {
            if (run != 0)
                pA.randomReset(A, model);
            pA.initialPartitioning(iters);

            std::vector<TIdx> communicationVolumes;
                communicationVolumes.push_back(A.communicationVolume());

            for (int i = 1; i <= 1000; ++i) {
                pA.refineWithIterations(A.nonZeros(), 0, 0, randomize);
                auto comVol = A.communicationVolume();
                communicationVolumes.push_back(comVol);
                if (communicationVolumes[i] == communicationVolumes[i - 1])
                    break;
            }

            ZeeLogVar(communicationVolumes);
            Vs.push_back((double)communicationVolumes.back());

            //A.spy("steam", true);
            if (plot) {
                auto p = Zee::Plotter<TIdx>();
                p["xlabel"] = "iterations";
                p["ylabel"] = "$V$";
                p["title"] = "Communication Volume PuLP ";
                p.addLine(communicationVolumes, "communicationVolumes");
                p.plot("comvol_pulp_" + matrix, true);
            }

        }

        double sum = std::accumulate(Vs.begin(), Vs.end(), 0.0, std::plus<double>());
        double average = sum / Vs.size();

        if (Vs.size() > 1) {
            double sd = 0.0;
            double squared_sum = std::accumulate(
                Vs.begin(), Vs.end(), 0.0, [average](double lhs, double rhs) {
                    return lhs + (rhs - average) * (rhs - average);
                });
            sd = sqrt(squared_sum / (Vs.size() - 1));
            report.addResult(matrix, "SD", sd);
        }

        report.addResult(matrix, "V_HP", average);

        ZeeLogVar(A.loadImbalance());

        CyclicPartitioner<decltype(A)> cyclic(procs);
        cyclic.partition(A);
        auto comVol = A.communicationVolume();
        improvements.push_back((comVol - average) / comVol);

        report.addResult(matrix, "V_C", comVol);

        if (compareMg && procs == 2) {
            MGPartitioner<decltype(A)> mg;
            mg.partition(A);
            while (!mg.locallyOptimal()) {
                ZeeLogVar(A.communicationVolume());
                mg.refine(A);
            }

            report.addResult(matrix, "V_MG", A.communicationVolume());
        }
    }

    if (!benchmark) {
        bench.silence();
    }
    bench.finish();

    report.print();

    double sumImprovements = std::accumulate(
        improvements.begin(), improvements.end(), 0.0, std::plus<double>());
    double averageImprovement = sumImprovements / improvements.size();
    ZeeLogResult << "Average improvement: " << std::fixed
                 << std::setprecision(2) << averageImprovement * 100.0 << "%"
                 << endLog;

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
