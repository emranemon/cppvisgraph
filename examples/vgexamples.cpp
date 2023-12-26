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

#include "cppvisgraph/vis_graph.hpp"

#include <iostream>
#include <vector>
#include <chrono>

typedef cppvisgraph::Point Point;

class vgexamples
{
public:
    static void example_vg1(void)
    {
        // Create a vector of polygons
        std::vector<std::vector<cppvisgraph::Point>> polygons;
        polygons.push_back({Point(0.0,1.0), Point(3.0,1.0), Point(1.5,4.0)});
        polygons.push_back({Point(4.0,4.0), Point(7.0,4.0), Point(5.5,8.0)});
        polygons.push_back({Point(5.0, 7.0), Point(7.0, 4.0), Point(3.5, 18.0)});
        polygons.push_back({Point(7.0, 7.0), Point(7.0, 14.0), Point(13.5, 18.0)});

        auto start_time = std::chrono::high_resolution_clock::now();
        cppvisgraph::VisGraph g;
        g.build(polygons);
        auto path = g.shortest_path(Point(1.5, 0.0), Point(4.0, 6.0));
        auto finish_time = std::chrono::high_resolution_clock::now();
        auto microsecs = std::chrono::duration_cast<std::chrono::microseconds>(finish_time - start_time);
        std::cout<<"Calc Time: "<<microsecs.count()<<" \u03BCs.\n";
        for(auto p: path)
        {
            std::cout<< p << " ";
        }
        std::cout<<std::endl;
    }
    
    static void example_vg2(void)
    {
        // Create a vector of polygons
        std::vector<std::vector<cppvisgraph::Point>> polygons;
        polygons.push_back({Point(1.148, 3.259), Point(1.540, 3.456), Point(2.267, 3.352), Point(2.546, 2.608), Point(4.299, -0.154), Point(3.892, -0.823), Point(3.389, -1.515), Point(2.081, -1.595), Point(1.897, -1.633), Point(1.179, 0.787), Point(0.690, 1.392), Point(0.690, 1.598), Point(0.990, 3.102), Point(1.148, 3.259)});
        polygons.push_back({Point(-2.738, 0.775), Point(-2.738, 1.458), Point(-2.660, 1.850), Point(-2.509, 2.303), Point(-2.271, 2.778), Point(-1.719, 3.331), Point(-1.719, 2.540), Point(-1.790, 1.688), Point(-1.861, 0.903), Point(-2.046, -0.022), Point(-2.525, -0.022), Point(-2.660, 0.383), Point(-2.738, 0.775)});

        auto start_time = std::chrono::high_resolution_clock::now();
        cppvisgraph::VisGraph g;
        g.build(polygons);
        auto path = g.shortest_path(Point(0.430, 1.068), Point(7.767, -9.629));
        auto finish_time = std::chrono::high_resolution_clock::now();
        auto microsecs = std::chrono::duration_cast<std::chrono::microseconds>(finish_time - start_time);
        std::cout<<"Calc Time: "<<microsecs.count()<<" \u03BCs.\n";
        for(auto p: path)
        {
            std::cout<< p << " ";
        }
        std::cout<<std::endl;
    }
};

int main()
{
    vgexamples::example_vg2();
    return 0;
}