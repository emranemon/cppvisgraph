// The MIT License (MIT)

// Copyright (c) 2023 MIEA MD EMON 
// https://github.com/emranemon

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "cppvisgraph/vis_graph.hpp"

namespace py = pybind11;

class Point2D: public cppvisgraph::Point 
{
public:
    using Point::Point;

    double get_x() const { return x(); }
    void set_x(){}
    double get_y() const { return y(); }
    void set_y(){}
};

class PyVisGraph 
{
private:
    cppvisgraph::VisGraph g;
public:
    void build(const std::vector<std::vector<Point2D>>& input)
    {
        const auto& points = reinterpret_cast<const std::vector<std::vector<cppvisgraph::Point>>&>(input);
        g.build(points);
    }
    
    std::vector<Point2D> shortest_path(const Point2D& origin, const Point2D& destination)
    {
        const auto& orgn = reinterpret_cast<const cppvisgraph::Point&>(origin);
        const auto& dest = reinterpret_cast<const cppvisgraph::Point&>(destination);
        std::vector<Point2D> path2D;
        auto path = g.shortest_path(orgn, dest);
        for (auto point : path)
        {
            Point2D point2D(point.x(), point.y());
            path2D.push_back(point2D);
        }
        return path2D;
    }
};

PYBIND11_MODULE(cppyvisgraph, m)
{
    py::class_<Point2D>(m, "Point")
        .def(py::init<double, double>())
        .def_property("x", &Point2D::get_x, &Point2D::set_x)
        .def_property("y", &Point2D::get_y, &Point2D::set_y)
        ;

    py::class_<PyVisGraph>(m, "VisGraph")
        .def(py::init<>())
        .def("build", &PyVisGraph::build)
        .def("shortest_path", &PyVisGraph::shortest_path)
        ;
}
