#include <map>
#include <string>
#include <sstream>
#include <iomanip>
#include <fstream>

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
        ss << std::fixed << std::setprecision(2) << result;
        entries_[row][column] = ss.str();

        if (ss.str().size() > columnWidth_[column]) {
            columnWidth_[column] = ss.str().size();
        }
    }

    void print() {
        ZeeLogResult << title_ << endLog;

        unsigned int lineSize = rowSize_ + 4;
        for (auto col : columnWidth_) {
            lineSize += col.second + 2;
        }
        std::string hline = "";
        for (unsigned int i = 0; i < lineSize; ++i)
            hline.push_back('-');

        auto addElement = [](int width, std::stringstream& result,
                             std::string entry) {
            result << std::left << std::setprecision(2) << std::setw(width) << std::setfill(' ')
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

    void saveToTex(std::string filename) {
        std::ofstream fout(filename);
        fout << "\\begin{table}" << std::endl;
        fout << "\\begin{tabular}{|l|";
        for (unsigned int i = 0; i < columns_.size(); ++i) {
            fout << "c";
            fout << ((i < (columns_.size() - 1)) ? " " : "|}");
        }
        fout << std::endl << "\\hline" << std::endl;
        fout << rowTitle_ << " & ";
        for (unsigned int i = 0; i < columns_.size();  ++i) {
            fout << "\\verb|" << columns_[i] << "|";
            fout << ((i < (columns_.size() - 1)) ? " & " : "\\\\");
        }
        fout << std::endl;
        fout << "\\hline" << std::endl;

        for (auto& rowCols : entries_) {
            fout << "\\verb|" << rowCols.first << "| & ";
            for (unsigned int i = 0; i < columns_.size();  ++i) {
                fout << rowCols.second[columns_[i]];
                fout << ((i < (columns_.size() - 1)) ? " & " : "\\\\");
            }
            fout << std::endl;
        }
        fout << "\\hline" << std::endl;
        fout << "\\end{tabular}" << std::endl;;

        fout << "\\caption{\\ldots}" << std::endl;
        fout << "\\end{table}" << std::endl;
    }

  private:
    std::string title_;
    std::string rowTitle_;
    std::map<std::string, std::map<std::string, std::string>> entries_;
    std::vector<std::string> columns_;
    std::map<std::string, unsigned int> columnWidth_;
    unsigned int rowSize_;
};
