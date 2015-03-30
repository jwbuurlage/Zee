/*
File: include/sparse_matrix.h

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.
*/

#pragma once

#include <string>

namespace Zee
{

using std::string;

enum class Color
{
    red,
    blue,
    yellow,
    green,
    clear
};

string colorOutput(Color tc) {
    switch(tc)
    {
        case Color::red:
            return "\033[1;31m";

        case Color::green:
            return "\033[1;32m";

        case Color::yellow:
            return "\033[1;33m";

        case Color::blue:
            return "\033[1;34m";

        case Color::clear:
        default:
            return "\033[0m";
    }
}

}
