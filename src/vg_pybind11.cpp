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
#include <iomanip>
#include <sstream>

#include "cppvisgraph/vis_graph.hpp"

namespace py = pybind11;

using namespace cppvisgraph;

std::string POINT_TO_STRING(const Point &p)
{
    std::ostringstream oss;
    oss << "Point(" << std::fixed << std::setprecision(2)
            << p.x() << ", " << p.y() << ")";
    return oss.str();
}

PYBIND11_MODULE(cppyvisgraph, m) 
{
    py::class_<Point>(m, "Point")
        .def(py::init<double, double>())
        .def_property_readonly("x", &Point::x)
        .def_property_readonly("y", &Point::y)
        // .def_property_readonly("polygon_id", &Point::polygon_id)
        .def("__repr__", [](const Point &p) {
            return POINT_TO_STRING(p);
        });   

    // py::class_<Edge>(m, "Edge")
    //     .def(py::init<const Point &, const Point &>())
    //     .def_property_readonly("p1", &Edge::p1)
    //     .def_property_readonly("p2", &Edge::p2)
    //     .def("__repr__", [](const Edge &e) {
    //         return "Edge(" + POINT_TO_STRING(e.p1()) + ", " + POINT_TO_STRING(e.p2()) + ")";
    //     });

    // py::class_<Graph>(m, "Graph")
    //     .def(py::init<>())
    //     .def(py::init<const std::vector<std::vector<Point>> &>())
    //     .def_property_readonly("edges", &Graph::get_edges)
    //     .def_property_readonly("polygons", &Graph::get_polygons);

    py::class_<VisGraph>(m, "VisGraph")
        .def(py::init<>())
        // .def("load", &VisGraph::load)
        // .def("save", &VisGraph::save)
        // .def("find_visible", &VisGraph::find_visible)
        // .def("point_in_polygon", &VisGraph::point_in_polygon)
        // .def("closest_point", &VisGraph::closest_point)
        // .def("update", &VisGraph::update, py::arg("points"), py::arg("origin") = nullptr, py::arg("destination") = nullptr)
        .def("build", &VisGraph::build, py::arg("input"), py::arg("multi_thread") = true, py::arg("status") = false)
        .def("shortest_path", &VisGraph::shortest_path);
}