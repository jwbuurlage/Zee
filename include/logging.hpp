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

#include <iostream>
#include <string>

#include "color_output.hpp"

namespace Zee {

using std::cerr;
using std::cout;
using std::endl;

void logError(string s) {
    cerr << colorOutput(Color::red) << 
        "ERROR: " << colorOutput(Color::clear) << s << endl;
}

void logInfo(string s) {
    cout << colorOutput(Color::yellow) << "INFO: " << 
        colorOutput(Color::clear) << s << endl;
}

class Logger {
    public:
        Logger();
        ~Logger();

        operator Logger <<(std::string rhs) {

        }
};

} // namespace Zee
