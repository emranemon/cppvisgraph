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

#ifndef CPPVISGRAPH_GRAPH_HPP
#define CPPVISGRAPH_GRAPH_HPP

#include <iostream>
#include <vector>
#include <limits>
#include <cmath>
#include <unordered_set>
#include <unordered_map>

namespace cppvisgraph 
{

class Point 
{
public:
    /* Point(double x, double y, int polygon_id = -1)
        : x_(std::round(x * 1e15) / 1e15), y_(std::round(y * 1e15) / 1e15), polygon_id_(polygon_id) {} */
    Point(double x, double y, int polygon_id = -1)
        : x_(x), y_(y), polygon_id_(polygon_id) {}
    
    /* bool operator==(const Point& point) const 
    {
        return std::fabs(x_ - point.x_) < std::numeric_limits<double>::epsilon()
            && std::fabs(y_ - point.y_) < std::numeric_limits<double>::epsilon();
    } */
    bool operator==(const Point& point) const 
    {
        return x_ == point.x_ && y_ == point.y_;
    }

    bool operator!=(const Point& point) const 
    {
        return !(*this == point);
    }

    double x() const { return x_; }

    double y() const { return y_; }

    int polygon_id() const { return polygon_id_; }

    void set_polygon_id(int pid) { polygon_id_ = pid; }

    friend std::ostream& operator<<(std::ostream& os, const Point& point)
    {
        os << "Point(" << point.x_ << ", " << point.y_ << ")";
        return os;
    }

private:
    double x_;
    double y_;
    int polygon_id_;
};

struct PointHash 
{
    std::size_t operator()(const Point& p) const 
    {
        return std::hash<double>()(p.x()) ^ std::hash<double>()(p.y());
    }
};

class Edge
{
public:
    Edge(const Point& point1, const Point& point2)
        : p1_(point1), p2_(point2) {}

    Point get_adjacent(const Point& point) const
    {
        if (point == p1_)
            return p2_;
        return p1_;
    }

    Point p1() const { return p1_; }

    Point p2() const { return p2_; }

    bool contains(const Point& point) const 
    {
        return p1_ == point || p2_ == point;
    }

    bool operator==(const Edge& edge) const 
    {
        return (p1_ == edge.p1_ && p2_ == edge.p2_) ||
               (p1_ == edge.p2_ && p2_ == edge.p1_);
    }

    bool operator!=(const Edge& edge) const 
    {
        return !(*this == edge);
    }

    friend std::ostream& operator<<(std::ostream& os, const Edge& edge)
    {
        os << "Edge(" << edge.p1_ << ", " << edge.p2_ << ")";
        return os;
    }
    
private:
    Point p1_;
    Point p2_;
};

struct EdgeHash 
{
    std::size_t operator()(const Edge& e) const 
    {
        return (std::hash<double>()(e.p1().x()) ^ std::hash<double>()(e.p1().y()))^
            (std::hash<double>()(e.p2().x()) ^ std::hash<double>()(e.p2().y()));
    }
};

class Graph 
{
public:
    Graph() = default;
    
    Graph(const std::vector<std::vector<Point>>& polygons)
    {
        int pid = 0;
        for (auto polygon : polygons)
        {
            if (polygon.front() == polygon.back() && polygon.size() > 1)
            {
                polygon.pop_back();
            }
            for (size_t i = 0; i < polygon.size(); ++i)
            {
                auto& point = polygon[i];
                auto& sibling_point = polygon[(i + 1) % polygon.size()];
                if (polygon.size() > 2)
                {
                    point.set_polygon_id(pid);
                    sibling_point.set_polygon_id(pid);
                    Edge edge(point, sibling_point);
                    polygons_[pid].insert(edge);
                    add_edge(edge);
                }
                else
                {
                    add_edge(Edge(point, sibling_point));
                }
            }
            if (polygon.size() > 2)
            {
                ++pid;
            }
        }
    }

    std::vector<Point> get_adjacent_points(const Point& point) const 
    {
        std::vector<Point> adjacent_points;
        if(graph_.count(point)>0)
        {
            const auto& edges = graph_.at(point);
            for (const auto& edge : edges)
            {
                adjacent_points.push_back(edge.get_adjacent(point));
            }
        }
        return adjacent_points;
    }

    std::vector<Point> get_points() const 
    {
        std::vector<Point> points;
        for (const auto& entry : graph_)
        {
            points.push_back(entry.first);
        }
        return points;
    }

    std::unordered_set<Edge, EdgeHash> get_edges() const 
    {
        return edges_;
    }

    std::unordered_map<int, std::unordered_set<Edge, EdgeHash>> get_polygons() const 
    {
        return polygons_;
    }

    std::unordered_set<Edge, EdgeHash> operator[](const Point& p) const 
    {
        std::unordered_set<Edge, EdgeHash> edges;
        if(graph_.count(p)>0)
            edges = graph_.at(p);
        return edges;
    }

    void add_edge(const Edge& edge)
    {
        graph_[edge.p1()].insert(edge);
        graph_[edge.p2()].insert(edge);
        edges_.insert(edge);
    }

    bool contains(const Point& point) const 
    {
        return graph_.count(point) > 0;
    }

    bool contains(const Edge& edge) const 
    {
        return edges_.count(edge) > 0;
    }

private:
    std::unordered_map<Point, std::unordered_set<Edge, EdgeHash>, PointHash> graph_;
    std::unordered_set<Edge, EdgeHash> edges_;
    std::unordered_map<int, std::unordered_set<Edge, EdgeHash>> polygons_;
};

}  // namespace cppvisgraph
#endif //CPPVISGRAPH_GRAPH_HPP