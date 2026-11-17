// Regression test: VisGraph::build()'s visibility sweep must never admit an edge between
// two non-adjacent vertices of the SAME convex obstacle polygon as "visible" - such a chord
// always cuts through that polygon's own interior. Root cause was a ray-cast degeneracy in
// PolygonUtils::polygon_crossing() (fixed alongside this test) that silently misclassified a
// midpoint that's genuinely inside a polygon as outside, whenever the ray happened to pass
// near/through an unrelated polygon vertex - a case regular N-gons hit often, since many of
// their own vertices land exactly on axis-aligned coordinates.
//
// Builds a Graph directly (same way VisGraph::build() does internally) and queries
// VisibleVertices::visible_vertices() on its own properly polygon_id-tagged points -
// VisGraph::find_visible() takes a caller-supplied Point by value, and Graph's constructor
// only tags its OWN internal copies (`for (auto polygon : polygons)` copies each inner
// vector), never the caller's originals, so querying with external, never-tagged points
// would silently skip the very check this test exists to exercise.
//
// No test framework dependency - this project doesn't have one. Exits 0 on success (CTest's
// pass convention), non-zero with a diagnostic on failure.
#include "cppvisgraph/graph.hpp"
#include "cppvisgraph/visible_vertices.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using namespace cppvisgraph;

namespace {

std::vector<Point> make_regular_polygon(double center_x, double center_y, double radius,
                                        int num_sides) {
    std::vector<Point> polygon;
    for (int i = 0; i < num_sides; ++i) {
        double angle = 2.0 * M_PI * i / num_sides;
        polygon.emplace_back(center_x + radius * std::cos(angle), center_y + radius * std::sin(angle));
    }
    return polygon;
}

bool approx_eq(const Point &a, const Point &b) {
    return std::abs(a.x() - b.x()) < 1e-6 && std::abs(a.y() - b.y()) < 1e-6;
}

} // namespace

int main() {
    int failures = 0;

    // A single convex 16-gon obstacle - no other library code involved. The chord between
    // vertex[0] and any non-adjacent vertex (every vertex except its two immediate
    // neighbors) must never be reported visible.
    constexpr int kSides = 16;
    std::vector<Point> polygon = make_regular_polygon(0.0, 0.0, 20.0, kSides);

    std::vector<std::vector<Point>> obstacles = {polygon};
    Graph graph(obstacles);

    // Properly polygon_id-tagged points, read back from the Graph itself - not our own
    // untagged originals.
    std::vector<Point> tagged = polygon;
    for (const auto &p : graph.get_points()) {
        for (int i = 0; i < kSides; ++i) {
            if (approx_eq(p, polygon[i])) tagged[i] = p;
        }
    }

    bool found_any_adjacent_edge = false;
    for (int i = 0; i < kSides; ++i) {
        std::vector<Point> visible_from_i = VisibleVertices::visible_vertices(tagged[i], graph);
        for (int j = 0; j < kSides; ++j) {
            if (i == j) continue;
            bool is_visible = false;
            for (const auto &v : visible_from_i) {
                if (approx_eq(v, tagged[j])) is_visible = true;
            }
            int step = std::abs(i - j);
            bool adjacent = (step == 1) || (step == kSides - 1);
            if (adjacent) {
                found_any_adjacent_edge = found_any_adjacent_edge || is_visible;
                continue; // adjacent pairs are legitimately visible - not this test's concern
            }
            if (is_visible) {
                std::fprintf(stderr,
                            "FAIL: same-polygon diagonal reported visible: vertex[%d]=(%.6f,%.6f) "
                            "-- vertex[%d]=(%.6f,%.6f) - this chord cuts through the polygon's "
                            "own interior and must never be a visibility-graph edge\n",
                            i, polygon[i].x(), polygon[i].y(), j, polygon[j].x(), polygon[j].y());
                ++failures;
            }
        }
    }

    // Sanity check the test itself isn't vacuous: a fix that over-corrected to "nothing on
    // this polygon is ever visible" would pass the loop above trivially without actually
    // testing anything - at least one legitimate adjacent-vertex (boundary) edge must exist.
    if (!found_any_adjacent_edge) {
        std::fprintf(stderr, "FAIL: no adjacent-vertex (legitimate boundary) edge was found at "
                             "all - the fix over-corrected\n");
        ++failures;
    }

    if (failures == 0) {
        std::printf("PASS: all %d same-polygon diagonals correctly rejected; boundary edges intact\n",
                   kSides);
        return 0;
    }
    std::fprintf(stderr, "%d failure(s)\n", failures);
    return 1;
}
