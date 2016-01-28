/*
File: include/plotter.hpp

This file is part of the Zee partitioning framework

Copyright (C) 2015 Jan-Willem Buurlage <janwillembuurlage@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License (LGPL)
as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This file has been adapted from the Arya game engine.
*/

#pragma once

#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <map>

#include "logging.hpp"

namespace Zee {

template <typename T>
struct Line {
    std::vector<T> ys;
    std::string title;
    std::string descriptor;
};

template <typename T = double>
class Plotter {
    public:
        void plot(std::string plotName = "anonymous", bool show = false) {
            constexpr const char* extension = ".yaml";

            // write .zplot file
            auto depth = 0;

            using std::endl;

            std::stringstream ss;
            ss << "data/plots/" << plotName << extension;
            auto filename = ss.str();
            int i = 1;
            while(fileExists(filename)) {
                ss.str("");
                ss.clear();
                ss << "data/plots/" << plotName << "_" << i++ << extension;
                filename = ss.str();
            }
            std::ofstream fout(filename);

            fout << "# file description: zee plot" << endl;

            for (const auto& pKeyValue : attributes_) {
                fout << pKeyValue.first << ": " <<
                    literalize(pKeyValue.second) << endl;
            }

            if (!lines_.empty())
                fout << "lines:" << endl;

            depth++;

            for (const auto& line : lines_) {
                fout << tabs(depth) << line.descriptor << ":" << endl;
                depth++;
                fout << tabs(depth) << "title: " << literalize(line.title) << endl;
                fout << tabs(depth) << "data:" << endl;
                depth++;
                for (const auto& val : line.ys) {
                    fout << tabs(depth) << "- " << val << endl;
                }
                depth -= 2;
            }

            if (show) {
                showFile(filename);
            }
        }


        void addLine(vector<T> data, std::string title)
        {
            Line<T> l{};

            std::stringstream ss;
            ss << "line_" << lines_.size();

            l.descriptor = ss.str();
            l.ys = data;
            l.title = title;
            lines_.push_back(l);
        }

        void showFile(std::string filename) const {
            auto command = "./script/plot.py --save --showfile " + filename;
            std::system(command.c_str());
        }

        std::string& operator[] (std::string attr) {
            return attributes_[attr];
        }

        // TODO: implement
        void reset() {

        }

    private:
        std::string tabs(int depth) const
        {
            std::stringstream tabs;
            for (int i = 0; i < depth; ++i) {
                tabs << "    ";
            }
            return tabs.str();
        }

        std::string literalize(std::string str) const
        {
            str = "\"" + str + "\"";
            return str;
        }

        std::vector<Line<T>> lines_;
        std::map<std::string, std::string> attributes_;
};

} // namespace Zee
