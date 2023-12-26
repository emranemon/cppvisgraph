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

#ifndef CPPVISGRAPH_SHORTEST_PATH_HPP
#define CPPVISGRAPH_SHORTEST_PATH_HPP

#include "graph.hpp"
#include "visible_vertices.hpp"

#include <unordered_map>
#include <vector>
#include <algorithm>
#include <functional>
#include <utility>
#include <queue>
#include <string>

namespace cppvisgraph
{

/**
* @brief Hashmap that can be used as a priority queue.
* This is a replica of priority_dict(dict) implemented in python
*
* Keys of the hashmap are items to be put into the queue, and values
* are their respective priorities. All hashmap methods work as expected.
* The advantage over a standard heapq-based priority queue is that priorities
* of items can be efficiently updated (amortized O(1)) using code as
* 'themap.push(item, new_priority).'
*
* Note that this is a modified version of
* https://gist.github.com/matteodellamico/4451520
*/
template <typename Key, typename Value, typename Hash = std::hash<Key>>
class priority_dict : public std::unordered_map<Key, Value, Hash>
{
private:
    typedef std::pair<Value, Key> ValueKeyPair;
    struct Greater
    {
        bool operator()(const ValueKeyPair& lhs, const ValueKeyPair& rhs)
        {
            return lhs.first > rhs.first;  // Greater than for min-heap, Smaller than for max-heap
        }
    };
    std::priority_queue<ValueKeyPair, std::vector<ValueKeyPair>, Greater> stl_heap;

public:
    void push(Key key, Value value)
    {
        auto result = this->emplace(key, value);
        if(!result.second)
        {
            result.first->second = value;
        }
        stl_heap.push(std::make_pair(value, key));
    }

    std::pair<Key, Value> pop()
    {
        auto vk = stl_heap.top();
        auto v = vk.first;
        auto k = vk.second;
        stl_heap.pop();
        while (this->count(k)==0 || this->at(k) != v)
        {
            vk = stl_heap.top();
            v = vk.first;
            k = vk.second;
            stl_heap.pop();
        }        
        this->erase(k);
        return std::pair<Key, Value>(k, v);
    }
}; // class priority_dict


class ShortestPath
{
private:
    static std::pair<std::unordered_map<Point, double, PointHash>, std::unordered_map<Point, Point, PointHash>> dijkstra(
            const Graph& graph, const Point& origin, const Point& destination, const Graph& add_to_visgraph)
    {
        std::unordered_map<Point, double, PointHash> D;
        std::unordered_map<Point, Point, PointHash> P;
        priority_dict<Point, double, PointHash> Q;
        Q.push(origin, 0);
        while (!Q.empty())
        {
            auto top = Q.pop();
            Point v = top.first;
            D[v] = top.second;
            if (v == destination) break;            
            std::unordered_set<Edge, EdgeHash> edges = graph[v];
            if (add_to_visgraph.contains(v) && !add_to_visgraph[v].empty())
            {
                edges.insert(add_to_visgraph[v].begin(), add_to_visgraph[v].end());
            }
            for (const Edge& e : edges)
            {
                Point w = e.get_adjacent(v);
                double elength = D.at(v) + PolygonUtils::edge_distance(v, w);
                if (D.find(w) != D.end())
                {
                    if (elength < D.at(w))
                    {
                        throw std::runtime_error("Shorter path found, contradicting the assumption of shortest distances.");
                    }
                }
                else if (Q.find(w) == Q.end() || elength < Q.at(w))
                {
                    Q.push(w, elength);
                    auto p = P.emplace(w, v);
                    if (!p.second) { p.first->second = v; }
                }
            }
        }
        return std::make_pair(D, P); 
    }

public:
    static std::vector<Point> get(const Graph& graph, const Point& origin, const Point& destination, 
            const Graph& add_to_visgraph = Graph())
    {
        auto dp = ShortestPath::dijkstra(graph, origin, destination, add_to_visgraph);
        std::unordered_map<Point, double, PointHash> D = dp.first;
        std::unordered_map<Point, Point, PointHash> P = dp.second;
        std::vector<Point> path;
        Point back_to_origin = destination;
        while (true)
        {
            path.push_back(back_to_origin);
            if (back_to_origin == origin) break;
            back_to_origin = P.at(back_to_origin);
        }
        std::reverse(path.begin(), path.end());
        return path;
    }
}; // class ShortestPath

}  // namespace cppvisgraph
#endif //CPPVISGRAPH_SHORTEST_PATH_HPP
