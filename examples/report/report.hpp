#include <map>
#include <string>
#include <sstream>
#include <iomanip>

class Report {
  public:
    Report(std::string title, std::string rowTitle)
        : title_(title), rowTitle_(rowTitle) {
        rowSize_ = rowTitle_.size();
    }

    void addColumn(std::string col) {
        columns_.push_back(col);
        columnWidth_[col] = col.size();
    }

    void addRow(std::string row) {
        entries_[row] = std::map<std::string, std::string>();

        if (row.size() > rowSize_) {
            rowSize_ = row.size();
        }
    }

    template <typename T>
    void addResult(std::string row, std::string column, T result) {
        if (entries_.find(row) == entries_.end()) {
            ZeeLogError << "Trying to add result to non-existing row" << endLog;
            return;
        }
        std::stringstream ss;
        ss << result;
        entries_[row][column] = ss.str();

        if (ss.str().size() > columnWidth_[column]) {
            columnWidth_[column] = ss.str().size();
        }
    }

    void print() {
        ZeeLogResult << title_ << endLog;

        auto hline =
            "----------------------------------------------------------";

        auto addElement = [](int width, std::stringstream& result,
                             std::string entry) {
            result << std::left << std::setw(width) << std::setfill(' ')
                   << entry;
        };

        std::stringstream ss;
        addElement(rowSize_ + 2, ss, rowTitle_);
        ss << "| ";

        for (auto& col : columns_)
            addElement(columnWidth_[col] + 2, ss, col);

        ZeeLogInfo << hline << endLog;
        ZeeLogInfo << ss.str() << endLog;
        ZeeLogInfo << hline << endLog;

        for (auto& rowCols : entries_) {
            std::stringstream rowSs;
            addElement(rowSize_ + 2, rowSs, rowCols.first);
            rowSs << "| ";
            for (auto& col : columns_) {
                addElement(columnWidth_[col] + 2, rowSs, rowCols.second[col]);
            }
            ZeeLogInfo << rowSs.str() << endLog;
        }
        ZeeLogInfo << hline << endLog;
    }

    void saveToCSV();
    void readFromCSV();

    void saveToTex();

  private:
    std::string title_;
    std::string rowTitle_;
    std::map<std::string, std::map<std::string, std::string>> entries_;
    std::vector<std::string> columns_;
    std::map<std::string, unsigned int> columnWidth_;
    unsigned int rowSize_;
};