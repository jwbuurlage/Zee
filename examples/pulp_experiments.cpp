#include <zee.hpp>

#include <limits>
#include <map>
#include <string>

using namespace Zee;

int main(int argc, char* argv[]) {
    using TVal = float;
    using TIdx = uint64_t;

    // available models
    static std::map<HGModel, std::string> modelNames = {
        {HGModel::row_net, "row-net"},
        {HGModel::column_net, "column-net"},
        {HGModel::fine_grain, "fine-grain"}};

    // parse CLI args
    auto args = jw::ArgParse();
    args.addOptionWithDefault("--hg", "hypergraph model to use, should be "
                                      "chosen in set { row_net, column_net, "
                                      "fine_grain, auto }",
                              "auto");
    args.addOptionWithDefault("--procs", "number of parts", 2);
    args.addOptionWithDefault("--eps", "tolerance level for load imbalance",
                              0.03);
    args.addOptionWithDefault("--iters",
                              "number of iterations in initial phase", 1);
    args.addOptionWithDefault("--runs",
                              "number of partitioning runs to average over", 1);
    args.addOption("--matrices", "list of matrices to partition", true);
    args.addOption("--plot", "show plot of convergence");
    args.addOption("--spy", "save spy plots of best partitionings");
    args.addOption("--benchmark", "benchmark partitioning");
    args.addOption("--randomize", "visit vertices in random order");
    args.addOption("--show-eps", "visit vertices in random order");
    args.addOption("--compare-mg", "also apply medium-grain to matrices");
    args.addOption("--save-to-tex-table",
                   "file to write result as tex table to ");
    args.addOption("--fraction-plot",
                   "plot the fraction of matrices with a minimum gain");
    if (!args.parse(argc, argv))
        return -1;

    auto matrices = args.asList("--matrices");
    JWLogVar(matrices);

    HGModel model = HGModel::fine_grain;
    auto modelName = args.asString("--hg");
    bool autoModel = false;

    JWLogVar(modelName);

    if (modelName == "fine_grain") {
        model = HGModel::fine_grain;
    } else if (modelName == "row_net") {
        model = HGModel::row_net;
    } else if (modelName == "column_net") {
        model = HGModel::column_net;
    } else if (modelName == "auto") {
        autoModel = true;
    }

    TIdx procs = args.as<TIdx>("--procs");
    auto epsilon = args.as<double>("--eps");
    bool plot = args.wasPassed("--plot");
    bool showEps = args.wasPassed("--show-eps");
    bool benchmark = args.wasPassed("--benchmark");
    bool randomize = args.wasPassed("--randomize");
    bool compareMg = args.wasPassed("--compare-mg");
    bool spyMatrix = args.wasPassed("--spy");
    bool fractionPlot = args.wasPassed("--fraction-plot");
    TIdx runs = args.as<TIdx>("--runs");
    TIdx iters = args.as<TIdx>("--iters");

    // initialize report
    auto report = Report("Hyper-PuLP", "matrix");
    report.addColumn("m");
    report.addColumn("n");
    report.addColumn("N");
    report.addColumn("V_HP", "V_{\\text{HP}}");
    report.addColumn("V_HP_m", "V_{\\text{HP}}^{\\text{min}}");
    report.addColumn("V_C");
    report.addColumn("G");
    if (compareMg && procs == 2)
        report.addColumn("V_MG");
    if (autoModel)
        report.addColumn("model", "\\text{model}");
    if (showEps)
        report.addColumn("eps", "$\\epsilon$");

    auto bench = Zee::Benchmark("PuLP");

    // perform experiments
    std::vector<double> improvements;
    std::vector<double> gains;
    for (auto matrix : matrices) {
        report.addRow(matrix);
        bench.phase(matrix);

        // initialize the matrix from file
        DSparseMatrix<TVal, TIdx> B{"data/matrices/" + matrix + ".mtx", 1};

        auto cyclicComVol = 0;
        if (autoModel) {
            CyclicPartitioner<decltype(B)> cyclicRow(procs, CyclicType::row);
            cyclicRow.partition(B);
            auto comVolRow = B.communicationVolume();

            CyclicPartitioner<decltype(B)> cyclicCol(procs, CyclicType::column);
            cyclicCol.partition(B);
            auto comVolColumn = B.communicationVolume();

            if (comVolRow > comVolColumn) {
                // Better to keep columns together
                model = HGModel::row_net;
                cyclicComVol = comVolColumn;
            } else {
                model = HGModel::column_net;
                cyclicComVol = comVolRow;
            }
        } else {
            if (model == HGModel::row_net) {
                CyclicPartitioner<decltype(B)> cyclic(procs,
                                                      CyclicType::column);
                cyclic.partition(B);
                cyclicComVol = B.communicationVolume();
            } else {
                // Also for fine_grain and others
                CyclicPartitioner<decltype(B)> cyclic(procs, CyclicType::row);
                cyclic.partition(B);
                cyclicComVol = B.communicationVolume();
            }
        }

        JWLogVar(modelNames[model]);

        DSparseMatrix<TVal, TIdx> A{"data/matrices/" + matrix + ".mtx", 1};

        report.addResult(matrix, "m", A.getRows());
        report.addResult(matrix, "n", A.getCols());
        report.addResult(matrix, "N", A.nonZeros());

        auto v = DVector<TVal, TIdx>{A.getCols(), 1.0};
        auto u = DVector<TVal, TIdx>{A.getRows()};

        PulpPartitioner<decltype(A)> pA(A, procs, epsilon);
        pA.initialize(A, model);

        std::vector<double> Vs;
        double bestV = std::numeric_limits<double>::max();
        double bestEps = -1.0;

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

            JWLogVar(communicationVolumes);
            Vs.push_back((double)communicationVolumes.back());
            if (Vs.back() < bestV) {
                bestV = Vs.back();
                bestEps = A.loadImbalance();
                if (spyMatrix) {
                    A.spy(matrix + "_hyperpulp");
                }
            }

            if (plot) {
                auto p = Zee::Plotter<TIdx>();
                p["xlabel"] = "iterations";
                p["ylabel"] = "$V$";
                p["title"] = "Communication Volume PuLP ";
                p.addLine(communicationVolumes, "communicationVolumes");
                p.plot("comvol_pulp_" + matrix, true);
            }
        }

        double sum =
            std::accumulate(Vs.begin(), Vs.end(), 0.0, std::plus<double>());
        double average = sum / Vs.size();

        if (Vs.size() > 1) {
            double sd = 0.0;
            double squared_sum = std::accumulate(
                Vs.begin(), Vs.end(), 0.0, [average](double lhs, double rhs) {
                    return lhs + (rhs - average) * (rhs - average);
                });
            sd = sqrt(squared_sum / (Vs.size() - 1));

            std::stringstream volumeResult;
            volumeResult << std::fixed << std::setprecision(1) << average
                         << " +- " << sd;

            report.addResult(matrix, "V_HP", volumeResult.str());
        } else {
            report.addResult(matrix, "V_HP", average);
        }
        report.addResult(matrix, "V_HP_m", bestV);

        if (showEps) {
            std::stringstream epsResult;
            epsResult << std::fixed << std::setprecision(4) << bestEps - 1.0;
            report.addResult(matrix, "eps", epsResult.str());
        }

        improvements.push_back((cyclicComVol - average) / cyclicComVol);

        report.addResult(matrix, "V_C", cyclicComVol);

        auto G = (1.0 - ((double)average / cyclicComVol)) * 100.0;
        std::stringstream gainResult;
        gainResult << std::fixed << std::setprecision(1) << G << "%";

        gains.push_back(G);
        report.addResult(matrix, "G", gainResult.str());

        if (autoModel) {
            report.addResult(matrix, "model", modelNames[model],
                             "\\text{" + modelNames[model] + "}");
        }

        if (compareMg && procs == 2) {
            MGPartitioner<decltype(A)> mg;
            mg.partition(A);
            while (!mg.locallyOptimal()) {
                JWLogVar(A.communicationVolume());
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
    if (args.wasPassed("--save-to-tex-table"))
        report.saveToTex(args.asString("--save-to-tex-table"));

    double sumImprovements = std::accumulate(
        improvements.begin(), improvements.end(), 0.0, std::plus<double>());
    double averageImprovement = sumImprovements / improvements.size();
    JWLogResult << "Average improvement: " << std::fixed
                 << std::setprecision(2) << averageImprovement * 100.0 << "%"
                 << endLog;

    // optionally plot fractions of matrices with a minimal gain
    if (fractionPlot) {
        TIdx totalMatrices = matrices.size();
        std::sort(gains.begin(), gains.end());
        std::vector<double> fractions;
        TIdx k = 0;
        for (double i = 0; i < 100.0; i += 1.0) {
            while (k < gains.size() && gains[k] < i) {
                ++k;
            }
            fractions.push_back(((double)totalMatrices - k) /
                                (double)totalMatrices);
        }

        JWLogVar(fractions);

        auto p = Zee::Plotter<double>();
        p["xlabel"] = "$G$";
        p["ylabel"] = "fraction";
        p.addLine(fractions, "fractions");
        p.plot("fractions", true);
    }

    return 0;
}
