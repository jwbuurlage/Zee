#include <map>
#include <vector>
#include <string>
#include <sstream>

class ArgParse {
    public:
        ArgParse(int argc, char* argv[]) {
            std::string currentFlag;
            for (int i = 1; i < argc; ++i) {
                if (argv[i][0] == '-') {
                    currentFlag = std::string(argv[i]);
                    args_[std::string(argv[i])] = "";
                } else {
                    if (currentFlag.empty()) {
                        ZeeLogError << "ArgParse: argument before flag." << endLog;
                        return;
                    }

                    auto val = std::string(argv[i]);

                    if (args_[currentFlag].empty()) {
                        args_[currentFlag] = val;
                    } else {
                        args_[currentFlag] += " " + val;
                    }
                }
            }

            ZeeLogVar(args_);
        }

        bool wasPassed(std::string key) {
            if (args_.find(key) != args_.end()) {
                return true;
            }
            return false;
        }

        std::vector<std::string> asList(std::string key) {
            std::vector<std::string> result;
            if (!wasPassed(key))
                return result;

            std::stringstream ss(args_[key]);
            while (ss.good()) {
                std::string element;
                ss >> element;
                result.push_back(element);
            }

            return result;
        }

        std::string asString(std::string key) {
            if (!wasPassed(key))
                return "";

            return args_[key];
        }

        template <typename T>
        T as(std::string key) {
            if (!wasPassed(key))
                return 0;

            T result = 0;
            std::stringstream ss(args_[key]);
            while (ss.good()) {
                ss >> result;
            }

            return result;
        }

    private:
        std::map<std::string, std::string> args_;
};
