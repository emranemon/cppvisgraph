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

#ifndef CPPVISGRAPH_VISIBLE_VERTICES_HPP
#define CPPVISGRAPH_VISIBLE_VERTICES_HPP

#include "graph.hpp"

#include <cmath>
#include <vector>
#include <limits>
#include <algorithm>
#include <string>
#include <tuple>

#define INF 10000
#define COLLINEAR 0
#define CCW 1
#define CW -1
#define DOUBLE_MAX std::numeric_limits<double>::max()
#define DOUBLE_MIN std::numeric_limits<double>::min()
#define COLIN_TOLERANCE 10
#define T pow(10, COLIN_TOLERANCE)
#define T2 pow(10.0, COLIN_TOLERANCE)

namespace cppvisgraph
{

class PolygonUtils
{
public:
    /**
    * @brief Returns true if Point p1 is internal to the polygon.
    *
    * The polygon is defined by the Edges in poly_edges. The method uses the crossings
    * algorithm and takes into account edges that are collinear to p1.
    *
    * @param p1 The Point to check for internal position within the polygon.
    * @param poly_edges The unordered set of Edges defining the polygon.
    *
    * @return True if Point p1 is internal to the polygon; otherwise, false.
    */
    static bool polygon_crossing(const Point& p1, const std::unordered_set<Edge, EdgeHash>& poly_edges)
    {
        Point p2(INF, p1.y());
        int intersect_count = 0;
        for (const Edge& edge : poly_edges)
        {
            if (p1.y() < edge.p1().y() && p1.y() < edge.p2().y()) continue;
            if (p1.y() > edge.p1().y() && p1.y() > edge.p2().y()) continue;
            if (p1.x() > edge.p1().x() && p1.x() > edge.p2().x()) continue;
            // Deal with points collinear to p1
            bool edge_p1_collinear = (ccw(p1, edge.p1(), p2) == COLLINEAR);
            bool edge_p2_collinear = (ccw(p1, edge.p2(), p2) == COLLINEAR);
            if (edge_p1_collinear && edge_p2_collinear) continue;
            if (edge_p1_collinear || edge_p2_collinear)
            {
                const Point& collinear_point = edge_p1_collinear ? edge.p1() : edge.p2();
                if (edge.get_adjacent(collinear_point).y() > p1.y())
                {
                    intersect_count += 1;
                }
            } 
            else if (edge_intersect(p1, p2, edge))
            {
                intersect_count += 1;
            }
        }
        return (intersect_count % 2 != 0);
    }

    /**
    * @brief Returns true if the edge from Point p1 to Point p2 is interior to any polygon in the graph.
    *
    * @param p1 The starting Point of the edge.
    * @param p2 The ending Point of the edge.
    * @param graph The Graph containing polygons to check for edge interiority.
    *
    * @return True if the edge from Point p1 to Point p2 is interior to any polygon in the graph; otherwise, false.
    */
    static bool edge_in_polygon(const Point& p1, const Point& p2, const Graph& graph)
    {
        if (p1.polygon_id() != p2.polygon_id()) return false;
        if (p1.polygon_id() == -1 || p2.polygon_id() == -1) return false;
        Point mid_point((p1.x() + p2.x()) / 2, (p1.y() + p2.y()) / 2);
        return polygon_crossing(mid_point, graph.get_polygons().at(p1.polygon_id()));
    }

    /**
    * @brief Returns Int (Greater than -1) if the Point p is interior to any polygon in the graph.
    *
    * @param p The Point to check for interiority within any polygon.
    * @param graph The Graph containing polygons to check for point interiority.
    *
    * @return Greater than -1 if the Point p is interior to any polygon in the graph; otherwise, -1.
    */
    static int point_in_polygon(const Point& p, const Graph& graph)
    {
        for (int i = 0; i < graph.get_polygons().size(); ++i)
        {
            if (polygon_crossing(p, graph.get_polygons().at(i)))
            {
                return i;
            }
        }
        return -1;
    }

    /**
    * @brief Calculates and returns the unit vector from Point c to Point p.
    *
    * The unit vector is a vector with a length of 1 and points in the direction from Point c to Point p.
    *
    * @param c The starting Point of the vector.
    * @param p The ending Point of the vector.
    *
    * @return The unit vector from Point c to Point p.
    */
    static Point unit_vector(const Point& c, const Point& p)
    {
        double magnitude = edge_distance(c, p);
        return Point((p.x() - c.x()) / magnitude, (p.y() - c.y()) / magnitude);
    }

    /**
    * @brief Assumes Point p is interior to the polygon with the specified polygon_id.
    *        Returns the closest point c outside the polygon to p, where the distance
    *        from c to the intersect point from p to the edge of the polygon is length.
    *
    * @param p The Point assumed to be interior to the polygon.
    * @param graph The Graph containing the polygon.
    * @param polygon_id The ID of the polygon to which Point p is assumed to be interior.
    * @param length The distance from the closest point c to the intersect point on the polygon edge (default is 0.001).
    *
    * @return The closest Point c outside the polygon to Point p.
    */
    static Point closest_point(const Point& p, const Graph& graph, int polygon_id, double length = 0.001)
    {
        auto polygon_edges = graph.get_polygons().at(polygon_id);
        Point close_point(0.0, 0.0);
        Edge close_edge(close_point, close_point);
        double close_dist = DOUBLE_MAX;
        for (const Edge& e : polygon_edges) 
        {
            double num = ((p.x() - e.p1().x()) * (e.p2().x() - e.p1().x()) + (p.y() - e.p1().y()) * (e.p2().y() - e.p1().y()));
            double denom = pow(edge_distance(e.p1(), e.p2()), 2);
            double u = num / denom;
            Point pu(
                e.p1().x() + u * (e.p2().x() - e.p1().x()),
                e.p1().y() + u * (e.p2().y() - e.p1().y())
            );
            if (u < 0.0)
            {
                pu = e.p1();
            } 
            else if(u > 1.0 )
            {
                pu = e.p2();
            }
            double d = edge_distance(p, pu);
            if(d < close_dist)
            {
                close_dist = d;
                close_point = pu;
                close_edge = e;
            }
        }
        if(close_edge.contains(close_point))
        {
            Point c = close_edge.p1() == close_point ? close_edge.p1() : close_edge.p2();
            std::unordered_set<Edge, EdgeHash> edgeSet = graph[c];
            std::vector<Edge> edges(edgeSet.begin(), edgeSet.end());
            Point v1 = unit_vector(c, edges[0].get_adjacent(c));
            Point v2 = unit_vector(c, edges[1].get_adjacent(c));
            Point vsum = unit_vector(Point(0, 0), Point(v1.x() + v2.x(), v1.y() + v2.y()));
            Point close1(c.x() + vsum.x() * length, c.y() + vsum.y() * length);
            Point close2(c.x() - vsum.x() * length, c.y() - vsum.y() * length);
            if (point_in_polygon(close1, graph) == -1)
            {
                return close1;
            }
            else
            {
                return close2;
            }
        }
        else
        {
            Point v = unit_vector(p, close_point);
            return Point(close_point.x() + v.x()*length, close_point.y() + v.y()*length);
        }
    }

    /**
    * @brief Returns the Euclidean distance between two Points.
    *
    * @param p1 The first Point.
    * @param p2 The second Point.
    *
    * @return The Euclidean distance between Point p1 and Point p2.
    */
    static double edge_distance(const Point& p1, const Point& p2)
    {
        return sqrt(pow((p2.x() - p1.x()), 2) + pow((p2.y() - p1.y()), 2));
    }
    
    /**
    * @brief Returns the intersect Point where the edge from Point p1 to Point p2 intersects with the specified Edge.
    *
    * @param p1 The starting Point of the edge.
    * @param p2 The ending Point of the edge.
    * @param edge The Edge to check for intersection.
    *
    * @return The intersect Point where the edge from Point p1 to Point p2 intersects with the specified Edge,
    *         or nullptr if there is no intersection.
    */
    static Point* intersect_point(const Point& p1, const Point& p2, const Edge& edge)
    {
        if (edge.contains(p1)) return new Point(p1.x(), p1.y());
        if (edge.contains(p2)) return new Point(p2.x(), p2.y());        
        if (edge.p1().x() == edge.p2().x())
        {
            if (p1.x() == p2.x()) 
                return nullptr;
            double pslope = (p1.y() - p2.y()) / (p1.x() - p2.x());
            double intersect_x = edge.p1().x();
            double intersect_y = pslope * (intersect_x - p1.x()) + p1.y();
            return new Point(intersect_x, intersect_y);
        }        
        if (p1.x() == p2.x())
        {
            double eslope = (edge.p1().y() - edge.p2().y()) / (edge.p1().x() - edge.p2().x());
            double intersect_x = p1.x();
            double intersect_y = eslope * (intersect_x - edge.p1().x()) + edge.p1().y();
            return new Point(intersect_x, intersect_y);
        }        
        double pslope = (p1.y() - p2.y()) / (p1.x() - p2.x());
        double eslope = (edge.p1().y() - edge.p2().y()) / (edge.p1().x() - edge.p2().x());
        if (eslope == pslope) return nullptr;        
        double intersect_x = (eslope * edge.p1().x() - pslope * p1.x() + p1.y() - edge.p1().y()) / (eslope - pslope);
        double intersect_y = eslope * (intersect_x - edge.p1().x()) + edge.p1().y();
        return new Point(intersect_x, intersect_y);
    }

    /**
    * @brief Returns the Euclidean distance from Point p1 to the intersect point with the specified Edge.
    *
    * Assumes the line going from Point p1 to Point p2 intersects the Edge before reaching Point p2.
    *
    * @param p1 The starting Point of the line.
    * @param p2 The ending Point of the line.
    * @param edge The Edge to check for intersection.
    *
    * @return The Euclidean distance from Point p1 to the intersect point with the specified Edge.
    */
    static double point_edge_distance(const Point& p1, const Point& p2, const Edge& edge)
    {
        Point* ip = intersect_point(p1, p2, edge);
        double p_edge_dist = 0.0;
        if (ip)
        {
            p_edge_dist = edge_distance(p1, *ip);
        }
        delete ip;
        return p_edge_dist;
    }

    /**
    * @brief Returns the angle (in radians) of Point point from the center of the radial circle.
    *
    * The angle is calculated with respect to the positive x-axis, where the center is the origin.
    *
    *         ------p
    *         |    /
    *         |   /
    *         |  /
    *        c|a/
    *
    * @param center The center Point of the radial circle.
    * @param point The Point for which the angle is calculated.
    *
    * @return The angle (in radians) of Point point from the center of the radial circle.
    */
    static double angle(const Point& center, const Point& point)
    {
        double dx = point.x() - center.x();
        double dy = point.y() - center.y();
        if (dx == 0)
        {
            if (dy < 0) 
                return M_PI * 3 / 2;
            return M_PI / 2;
        }
        if (dy == 0)
        {
            if (dx < 0) 
                return M_PI;
            return 0;
        }
        if (dx < 0) return M_PI + atan(dy / dx);
        if (dy < 0) return 2 * M_PI + atan(dy / dx);
        return atan(dy / dx);
    }

    /**
    * @brief Returns the angle B (in radians) between Point point_b and Point point_c.
    *
    * The angle is calculated within the context of a triangle with vertices at Point point_a, point_b, and point_c.
    * Angle B is the angle at vertex point_b.
    *
    *        c
    *      /  \
    *    /    B \
    *   a-------b
    *
    * @param point_a The first vertex of the triangle.
    * @param point_b The second vertex of the triangle where angle B is calculated.
    * @param point_c The third vertex of the triangle.
    *
    * @return The angle B (in radians) between Point point_b and Point point_c.
    */
    static double angle2(const Point& point_a, const Point& point_b, const Point& point_c)
    {
        double a = pow(point_c.x() - point_b.x(), 2) + pow(point_c.y() - point_b.y(), 2);
        double b = pow(point_c.x() - point_a.x(), 2) + pow(point_c.y() - point_a.y(), 2);
        double c = pow(point_b.x() - point_a.x(), 2) + pow(point_b.y() - point_a.y(), 2);
        double cos_value = (a + c - b) / (2 * sqrt(a) * sqrt(c));
        return acos(round(cos_value*T)/T2);
    }

    /**
    * @brief Returns true if the edge from Point A to Point B intersects with the specified Edge.
    *
    * The method uses the algorithm from http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
    *
    * @param p1 The starting Point of the first edge.
    * @param q1 The ending Point of the first edge.
    * @param edge The Edge to check for intersection with the first edge.
    *
    * @return True if the edge from Point A to Point B intersects with the specified Edge; otherwise, false.
    */
    static bool edge_intersect(const Point& p1, const Point& q1, const Edge& edge)
    {
        const Point& p2 = edge.p1();
        const Point& q2 = edge.p2();
        int o1 = ccw(p1, q1, p2);
        int o2 = ccw(p1, q1, q2);
        int o3 = ccw(p2, q2, p1);
        int o4 = ccw(p2, q2, q1);
        if (o1 != o2 && o3 != o4) return true;
        if (o1 == COLLINEAR && on_segment(p1, p2, q1)) return true;
        if (o2 == COLLINEAR && on_segment(p1, q2, q1)) return true;
        if (o3 == COLLINEAR && on_segment(p2, p1, q2)) return true;
        if (o4 == COLLINEAR && on_segment(p2, q1, q2)) return true;
        return false;
    }

    /**
    * @brief Given three collinear points p, q, r, the function checks if point q lies on line segment 'pr'.
    *
    * @param p The first collinear Point.
    * @param q The Point to check if it lies on the line segment.
    * @param r The second collinear Point.
    *
    * @return True if Point q lies on the line segment 'pr'; otherwise, false.
    */
    static bool on_segment(const Point& p, const Point& q, const Point& r)
    {
        if (q.x() <= std::max(p.x(), r.x()) && q.x() >= std::min(p.x(), r.x()) &&
            q.y() <= std::max(p.y(), r.y()) && q.y() >= std::min(p.y(), r.y()))
        {
            return true;
        }
        return false;
    }

    /**
    * @brief Returns 1 if the points A, B, C form a counter-clockwise orientation,
    *        -1 if they form a clockwise orientation, and 0 if they are collinear.
    *
    * @param A The first Point.
    * @param B The second Point.
    * @param C The third Point.
    *
    * @return 1 if the points A, B, C form a counter-clockwise orientation,
    *        -1 if they form a clockwise orientation, and 0 if they are collinear.
    */
    static int ccw(const Point& A, const Point& B, const Point& C)
    {
        double _area = (B.x() - A.x()) * (C.y() - A.y()) - (B.y() - A.y()) * (C.x() - A.x());
        double area = round(_area*T)/T2;
        if (area > 0.0)
            return CCW;
        if (area < 0.0)
            return CW;
        return 0;
    }
};

class OpenEdges
{
public:
    OpenEdges() {}

    void insert(const Point& p1, const Point& p2, const Edge& edge)
    {
        _open_edges.insert(_open_edges.begin() + _index(p1, p2, edge), edge);
    }

    void remove(const Point& p1, const Point& p2, const Edge& edge)
    {
        int index = _index(p1, p2, edge) - 1;
        if (_open_edges[index] == edge)
        {
            _open_edges.erase(_open_edges.begin() + index);
        }
    }

    Edge smallest() const
    {
        return _open_edges[0];
    }

    std::size_t size() const
    {
        return _open_edges.size();
    }

    Edge operator[](std::size_t index) const
    {
        return _open_edges[index];
    }

    std::vector<Edge> get_edges() const
    {
        return _open_edges;
    }

private:
    std::vector<Edge> _open_edges;

    /**
    * @brief Returns true if Edge edge1 is smaller than Edge edge2; otherwise, returns false.
    *
    * The comparison is based on the angles formed by the edges with respect to the given Points p1 and p2.
    *
    * @param p1 The reference Point for angle calculation related to edge1 and edge2.
    * @param p2 The second reference Point for angle calculation.
    * @param edge1 The first Edge to compare.
    * @param edge2 The second Edge to compare.
    *
    * @return True if Edge edge1 is smaller than Edge edge2; otherwise, false.
    */
    bool _less_than(const Point& p1, const Point& p2, const Edge& edge1, const Edge& edge2) const
    {
        if (edge1 == edge2) { return false; }
        if (!PolygonUtils::edge_intersect(p1, p2, edge2)) { return true; }
        double edge1_dist = PolygonUtils::point_edge_distance(p1, p2, edge1);
        double edge2_dist = PolygonUtils::point_edge_distance(p1, p2, edge2);
        if (edge1_dist > edge2_dist) 
        { 
            return false;
        }
        else if (edge1_dist < edge2_dist)
        { 
            return true;
        }
        else
        {
            Point same_point(0.0, 0.0);
            if (edge2.contains(edge1.p1())) 
            {
                same_point = edge1.p1();
            } 
            else
            {
                same_point = edge1.p2();
            }
            double angle_edge1 = PolygonUtils::angle2(p1, p2, edge1.get_adjacent(same_point));
            double angle_edge2 = PolygonUtils::angle2(p1, p2, edge2.get_adjacent(same_point));
            if (angle_edge1 < angle_edge2)
            {
                return true;
            }
            return false;
        }
    }

    std::size_t _index(const Point& p1, const Point& p2, const Edge& edge) const
    {
        int lo = 0;
        int hi = _open_edges.size();
        while (lo < hi)
        {
            int mid = (lo + hi) / 2;
            if (_less_than(p1, p2, edge, _open_edges[mid]))
            {
                hi = mid;
            } 
            else
            {
                lo = mid + 1;
            }
        }
        return lo;
    }
};

class PointComparison
{
private:
    const Point& point;

public:
    PointComparison(const Point& p) : point(p) {}

    bool operator()(const Point& p1, const Point& p2) const
    {
        double angle1 = PolygonUtils::angle(point, p1);
        double distance1 = PolygonUtils::edge_distance(point, p1);
        double angle2 = PolygonUtils::angle(point, p2);
        double distance2 = PolygonUtils::edge_distance(point, p2);
        //return angle1 < angle2 || (angle1 == angle2 && distance1 < distance2);
        return std::tie(angle1, distance1) < std::tie(angle2, distance2);
    }
};

class VisibleVertices
{
public:
    /**
    * @brief Returns a list of Points in the graph visible by a specified point.
    *
    * If origin and/or destination Points are given, these will also be checked for visibility.
    * The scan parameter can be set to 'full' to check for visibility against all points in the graph,
    * or 'half' to check for visibility against half the points. This saves running time when building
    * a complete visibility graph, as the points that are not checked will eventually be 'point.'
    *
    * @param point The reference Point from which visibility is checked.
    * @param graph The Graph containing Points and Edges to check visibility against.
    * @param origin Optional origin Point to check for visibility.
    * @param destination Optional destination Point to check for visibility.
    * @param scan Scan mode, either 'full' or 'half.'
    *
    * @return List of Points in the graph visible by the specified point.
    */
    static std::vector<Point> visible_vertices(const Point& point, const Graph& graph, const Point* origin = nullptr, const Point* destination = nullptr, const std::string& scan = "full")
    {
        auto edges = graph.get_edges();
        std::vector<Point> points = graph.get_points();
        if (origin) points.push_back(*origin);
        if (destination) points.push_back(*destination);        
        std::sort(points.begin(), points.end(), PointComparison(point));
        OpenEdges open_edges;
        Point point_inf(INF, point.y());
        for (const Edge& edge : edges)
        {
            if (edge.contains(point)) continue;
            if (PolygonUtils::edge_intersect(point, point_inf, edge))
            {
                if (PolygonUtils::on_segment(point, edge.p1(), point_inf)) continue;
                if (PolygonUtils::on_segment(point, edge.p2(), point_inf)) continue;
                open_edges.insert(point, point_inf, edge);
            }
        }
        std::vector<Point> visible;
        const Point* prev = nullptr;
        bool prev_visible = false;
        for (const Point& p : points)
        {
            if (p == point) continue;
            if (scan == "half" && PolygonUtils::angle(point, p) > M_PI) break;
            if (open_edges.size()!=0)
            {
                for (const Edge& edge : graph[p])
                {
                    if (PolygonUtils::ccw(point, p, edge.get_adjacent(p)) == CW)
                    {
                        open_edges.remove(point, p, edge);
                    }
                }
            }
            bool is_visible = false;
            if (prev == nullptr || PolygonUtils::ccw(point, *prev, p) != COLLINEAR || !PolygonUtils::on_segment(point, *prev, p))
            {
                if (open_edges.size()==0 || !PolygonUtils::edge_intersect(point, p, open_edges.smallest()))
                {
                    is_visible = true;
                }
            } 
            else if (!prev_visible)
            {
                is_visible = false;
            }
            else
            {
                is_visible = true;
                for (const Edge& edge : open_edges.get_edges())
                {
                    if (!edge.contains(*prev) && PolygonUtils::edge_intersect(*prev, p, edge))
                    {
                        is_visible = false;
                        break;
                    }
                }
                if (is_visible && PolygonUtils::edge_in_polygon(*prev, p, graph))
                {
                    is_visible = false;
                }
            }
            auto ap = graph.get_adjacent_points(point);            
            if (is_visible && std::find(ap.begin(), ap.end(), p) == ap.end())
            {
                is_visible = !PolygonUtils::edge_in_polygon(point, p, graph);
            }
            if (is_visible)
            {
                visible.push_back(p);
            }
            for (const Edge& edge : graph[p])
            {
                if (!edge.contains(point) && PolygonUtils::ccw(point, p, edge.get_adjacent(p)) == CCW)
                {
                    open_edges.insert(point, p, edge);
                }
            }
            prev = &p;
            prev_visible = is_visible;
        }
        return visible;
    }
};

}  // namespace cppvisgraph
#endif //CPPVISGRAPH_VISIBLE_VERTICES_HPP
