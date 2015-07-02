/*
File: include/logging.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This file has been adapted from the Arya game engine.
*/

#pragma once

#include <sstream>
#include <iostream>
#include <string>

#include "color_output.hpp"

#define ZeeInfoLog (Zee::Logger() << Zee::LogType::info)
#define ZeeWarningLog (Zee::Logger() << Zee::LogType::warning)
#define ZeeErrorLog (Zee::Logger() << Zee::LogType::error)

#define endLog Zee::Logger::end()

namespace Zee {

using std::cerr;
using std::cout;
using std::endl;

enum LogType {
    info,
    warning,
    error
};

class Logger {

    public:
        struct end { };

        Logger& operator <<(LogType t) {
            _t = t;
            return *this;
        }

        template <typename T>
        Logger& operator <<(T rhs) {
            ss << rhs;
            return *this;
        }

        void operator <<(end) {
            // output ss
            switch (_t) {
                case LogType::info:
                    cout << colorOutput(Color::yellow) << "INFO: ";
                    break;

                case LogType::warning:
                    cout << colorOutput(Color::blue) << "WARNING: ";
                    break;

                case LogType::error:
                    cout << colorOutput(Color::red) << "ERROR: ";
                    break;

                default:
                    break;
            }

            cout << colorOutput(Color::clear);
            cout << ss.str() << endl;
        }

    private:
        std::stringstream ss;
        LogType _t = LogType::info;
};


} // namespace Zee
