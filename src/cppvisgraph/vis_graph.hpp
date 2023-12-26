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

#ifndef CPPVISGRAPH_VIS_GRAPH_HPP
#define CPPVISGRAPH_VIS_GRAPH_HPP

#include "graph.hpp"
#include "visible_vertices.hpp"
#include "shortest_path.hpp"

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>
#include <queue>
#include <future>
#include <functional>

using namespace std::placeholders;

namespace cppvisgraph
{

class VisGraph
{
private:
    Graph graph;
    Graph visgraph;

public:
    void load(const std::string& filename)
    {
        // Load obstacle graph and visibility graph from file
        std::ifstream file(filename, std::ios::binary);
        if (file.is_open())
        {
            // Read the graph and visgraph from the file
            // Assuming the file format is compatible with the class structure
            file.read(reinterpret_cast<char*>(&graph), sizeof(Graph));
            file.read(reinterpret_cast<char*>(&visgraph), sizeof(Graph));
            file.close();
        } 
        else
        {
            // File open error
            std::cerr << "Error: Failed to open file " << filename << std::endl;
        }
    }

    void save(const std::string& filename)
    {
        // Save obstacle graph and visibility graph to file
        std::ofstream file(filename, std::ios::binary);
        if (file.is_open())
        {
            // Write the graph and visgraph to the file
            // Assuming the file format is compatible with the class structure
            file.write(reinterpret_cast<const char*>(&graph), sizeof(Graph));
            file.write(reinterpret_cast<const char*>(&visgraph), sizeof(Graph));
            file.close();
        } 
        else
        {
            // File open error
            std::cerr << "Error: Failed to open file " << filename << std::endl;
        }
    }

    /**
    * @brief Build visibility graph based on a list of polygons.
    *
    * This method constructs a visibility graph based on a list of polygons. The visibility graph
    * represents the visibility relationships between points in the given polygons.
    *
    * @param input List of polygons, where each polygon is represented by a vector of points.
    *              Each polygon should be a list of in-order (clockwise or counter-clockwise) Points.
    *              If there's only one polygon, it must still be a list within a list, e.g.,
    *              [[Point(0,0), Point(2,0), Point(2,1)]].
    * @param multi_thread Flag indicating whether to use multi-threading during the build process.
    * @param status Flag indicating whether to print status information during the build process.
    *
    * @note The implementation details of the visibility graph construction go here.
    */
    void build(const std::vector<std::vector<Point>>& input, const bool multi_thread = true, const bool status = false)
    {
        // Build visibility graph based on a list of polygons
        graph = Graph(input);
        visgraph = Graph();
        std::vector<Point> points = graph.get_points();
        const int batch_size = 10;
        if(status)
        {
            // Implement progress bar
            std::cout << "Progress bar is not implemented yet." << std::endl;
        }
        if (multi_thread)
        {
            std::vector<std::future<std::vector<Edge>>> futures;
            for (size_t i = 0; i < points.size(); i += batch_size)
            {
                futures.push_back(std::async(std::launch::async, std::bind(&VisGraph::_vis_graph, this, _1, _2), std::cref(graph),
                                  get_batch(points, i, batch_size)));
            }
            for (auto& future : futures)
            {
                for (const auto& edge : future.get())
                {
                    visgraph.add_edge(edge);
                }
            }
        } 
        else
        {
            for (size_t i = 0; i < points.size(); i += batch_size)
            {        
                std::vector<Point> batch = get_batch(points, i, batch_size);
                for (const auto& edge : _vis_graph(graph, batch))
                {
                    visgraph.add_edge(edge);
                }
            }
        }
    }

    std::vector<Point> find_visible(const Point& point)
    {
        // Find vertices visible from the given point
        return VisibleVertices::visible_vertices(point, graph);
    }

    void update(const std::vector<Point>& points, const Point* origin = nullptr, const Point* destination = nullptr)
    {
        // Update visgraph by checking visibility of points
        for (const auto& p : points)
        {
            for (const auto& v : VisibleVertices::visible_vertices(p, graph, origin, destination))
            {
                visgraph.add_edge(Edge(p, v));
            }
        }
    }

    /**
    * @brief Find and return the shortest path between origin and destination.
    *
    * This method calculates and returns the in-order list of Points representing
    * the shortest path between the given origin and destination Points within
    * the visibility graph.
    *
    * If the origin or destination Points are not part of the visibility graph,
    * their respective visibility edges will be found temporarily for path calculation.
    *
    * @param origin The starting Point of the path.
    * @param destination The destination Point of the path.
    *
    * @return In-order list of Points representing the shortest path found.
    */
    std::vector<Point> shortest_path(const Point& origin, const Point& destination)
    {
        bool origin_exists = visgraph.contains(origin);
        bool destination_exists = visgraph.contains(destination);
        if (origin_exists && destination_exists)
        {
            return ShortestPath::get(visgraph, origin, destination);
        }
        Point* orgn = origin_exists ? nullptr : new Point(origin.x(), origin.y());
        Point* dest = destination_exists ? nullptr : new Point(destination.x(), destination.y());
        Graph add_to_visg;
        if (!origin_exists)
        {
            for (const auto& v : VisibleVertices::visible_vertices(origin, graph, nullptr, dest))
            {                
                add_to_visg.add_edge(Edge(origin, v));
            }
        }
        if (!destination_exists)
        {
            for (const auto& v : VisibleVertices::visible_vertices(destination, graph, orgn))
            {
                add_to_visg.add_edge(Edge(destination, v));
            }
        }
        delete orgn;
        delete dest;
        return ShortestPath::get(visgraph, origin, destination, add_to_visg);
    }

    int point_in_polygon(const Point& point)
    {
        // Return the polygon_id if the point is inside a polygon, -1 otherwise
        return PolygonUtils::point_in_polygon(point, graph);
    }

    Point closest_point(const Point& point, int polygon_id, double length = 0.001)
    {
        // Return the closest point outside the polygon from the given point
        return PolygonUtils::closest_point(point, graph, polygon_id, length);
    }

private:
    std::vector<Edge> _vis_graph(const Graph& graph, const std::vector<Point>& points)
    {
        // Compute visible edges for a batch of points
        std::vector<Edge> visible_edges;
        for (const auto& p1 : points)
        {
            for (const auto& p2 : VisibleVertices::visible_vertices(p1, graph, nullptr, nullptr, "half"))
            {
                visible_edges.push_back(Edge(p1, p2));
            }
        }
        return visible_edges;
    }

    std::vector<Point> get_batch(const std::vector<Point>& points, const size_t index, const int batch_size)
    {
        // Slice points into batch
        return std::vector<Point>(points.begin() + index, points.begin() + std::min(index + batch_size, points.size()));
    }
};

}  // namespace cppvisgraph
#endif // CPPVISGRAPH_VIS_GRAPH_HPP