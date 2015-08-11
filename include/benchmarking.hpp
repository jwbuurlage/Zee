#pragma once

#include <chrono>
#include <vector>
#include <utility>
#include <sstream>

namespace Zee
{

class Benchmark
{
    public:
        using TTimePoint = std::chrono::time_point<std::chrono::high_resolution_clock>;

        Benchmark(std::string title)
            : title_(title)
        {
            start_ = std::chrono::high_resolution_clock::now();
        }

        ~Benchmark() {
            auto end = std::chrono::high_resolution_clock::now();

            std::stringstream splitOutput;
            if (!splits_.empty()) {
                splits_.push_back(make_pair("", end));
                auto hline = "--------------------------------------------------";
                splitOutput << "\n" << hline << '\n';
                for (unsigned int i = 0; i < splits_.size() - 1; ++i) {
                    auto splitTime = splits_[i + 1].second - splits_[i].second;
                    auto ms = std::chrono::duration<double, std::milli>(splitTime).count();
                    splitOutput << splits_[i].first << "\t" << ms << " ms" << endl;
                }
                splitOutput << hline;
            }

            ZeeLogBenchmark << title_ << " total runtime: " <<
                std::chrono::duration<double, std::milli>(end - start_).count() <<
                " ms" <<
                splitOutput.str() << endLog;
        }

        void phase(std::string splitTitle)
        {
            splitTitle.resize(30, ' ');
            auto now = std::chrono::high_resolution_clock::now();
            splits_.push_back(make_pair(splitTitle, now));
        }

    private:
        std::vector<std::pair<std::string, TTimePoint>> splits_;
        std::string title_;
        TTimePoint start_;
};

} // namespace Zee
