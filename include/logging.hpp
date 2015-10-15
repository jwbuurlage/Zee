/*
File: include/logging.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include <sstream>
#include <iostream>
#include <string>
#include <array>
#include <vector>
#include <set>
#include <map>

#include "color_output.hpp"

#define ZeeLogBenchmark (Zee::Logger() << Zee::LogType::benchmark)
#define ZeeLogDebug (Zee::Logger() << Zee::LogType::debug)
#define ZeeLogError (Zee::Logger() << Zee::LogType::error)
#define ZeeLogInfo (Zee::Logger() << Zee::LogType::info)
#define ZeeLogWarning (Zee::Logger() << Zee::LogType::warning)

#define endLog Zee::Logger::end()

#define ZeeLogVar(VAR) (Zee::Logger() << Zee::LogType::debug << #VAR " = " << VAR << endLog)

#ifndef DEBUG
#define ZeeAssert(ASSERT) if (!(ASSERT)) {\
    ZeeLogError << "assertion '" #ASSERT "' failed at " << __FILE__  << ":" << __LINE__ << endLog;\
    exit(-1);\
}
#define ZeeAssertMsg(ASSERT, MESSAGE) if (!(ASSERT)) {\
    ZeeLogError << "assertion '" #ASSERT "' failed at " << __FILE__  << ":" << __LINE__ << "\n" << "\"" << MESSAGE << "\"" << endLog;\
    exit(-1);\
}
#else
#define ZeeAssert(ASSERT)
#define ZeeAssertMsg(ASSERT, MESSAGE)
#endif

namespace Zee {

using std::cerr;
using std::cout;
using std::endl;

enum LogType {
    info,
    warning,
    error,
    debug,
    benchmark
};

class Logger {

    public:
        ~Logger() {
            if (!finalized_)
                finalize();
        }

        struct end { };

        Logger& operator <<(LogType t) {
            t_ = t;
            return *this;
        }

        template <typename T>
        Logger& operator <<(const T& rhs) {
            ss << rhs;
            return *this;
        }

        Logger& operator <<(const bool& rhs)
        {
            if (rhs)
                ss << "true";
            else
                ss << "false";
            return *this;
        }

        template <typename S>
        Logger& operator <<(const std::vector<S>& rhs) {
            auto sep = "";
            *this << "[";
            for (S value : rhs) {
                *this << sep << value;
                sep = ", ";
            }
            *this << "]";
            return *this;
        }

        template <typename S, size_t size>
        Logger& operator <<(const std::array<S, size>& rhs) {
            auto sep = "";
            *this << "[";
            for (S value : rhs) {
                *this << sep << value;
                sep = ", ";
            }
            *this << "]";
            return *this;
        }



        template <typename S>
        Logger& operator <<(const std::set<S>& rhs) {
            auto sep = "";
            *this << "[";
            for (S value : rhs) {
                *this << sep << value;
                sep = ", ";
            }
            *this << "]";
            return *this;
        }

        template <typename S, typename T>
        Logger& operator <<(const std::map<S, T>& rhs) {
            auto sep = "";
            *this << "{";
            for (auto& pair : rhs) {
                *this << sep << pair.first << " -> " << pair.second;
                sep = ", ";
            }
            *this << "}";
            return *this;
        }

        template <typename S, typename T>
        Logger& operator <<(std::pair<S, T> rhs) {
            ss << "[ " << rhs.first << ",\t" << rhs.second << " ]";
            return *this;
        }

        void operator <<(end) {
            finalize();
        }

    private:
        std::stringstream ss;
        LogType t_ = LogType::info;
        bool finalized_ = false;

        void finalize() {
            switch (t_) {
                case LogType::info:
                    cout << Zee::colors::start["cyan"] << "INFO: ";
                    break;

                case LogType::warning:
                    cout << Zee::colors::start["blue"] << "WARNING: ";
                    break;

                case LogType::error:
                    cout << Zee::colors::start["red"] << "ERROR: ";
                    break;

                case LogType::debug:
                    cout << Zee::colors::start["darkgray"] << "DEBUG: ";
                    break;

                case LogType::benchmark:
                    cout << Zee::colors::start["purple"] << "BENCHMARK: ";
                    break;

                default:
                    break;
            }

            cout << Zee::colors::end;
            cout << ss.str() << endl;

            finalized_ = true;
        }
};

} // namespace Zee
