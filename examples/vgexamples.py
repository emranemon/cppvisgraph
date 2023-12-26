# The MIT License (MIT)

# Copyright (c) 2023 MIEA MD EMON 
# https://github.com/emranemon

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import sys
import os
sys.path.insert(0, 
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import build.cppyvisgraph as vg
import time

def example_vg1():
    polys = [[vg.Point(0.0, 1.0), vg.Point(3.0, 1.0), vg.Point(1.5, 4.0)],
                [vg.Point(4.0, 4.0), vg.Point(7.0, 4.0), vg.Point(5.5, 8.0)],
                [vg.Point(5.0, 7.0), vg.Point(7.0, 4.0), vg.Point(3.5, 18.0)],
                [vg.Point(7.0, 7.0), vg.Point(7.0, 14.0), vg.Point(13.5, 18.0)]]
    start = time.perf_counter()
    g = vg.VisGraph()
    g.build(polys)
    path = g.shortest_path(vg.Point(1.5, 0.0), vg.Point(4.0, 6.0))
    elapsed_time = (time.perf_counter() - start) * 1e6
    print(f"Calc time: {elapsed_time:.2f} µs")
    print(path)

def example_vg2():
    polys = [[vg.Point(1.148, 3.259), vg.Point(1.540, 3.456), vg.Point(2.267, 3.352), vg.Point(2.546, 2.608), vg.Point(4.299, -0.154), vg.Point(3.892, -0.823), vg.Point(3.389, -1.515), vg.Point(2.081, -1.595), vg.Point(1.897, -1.633), vg.Point(1.179, 0.787), vg.Point(0.690, 1.392), vg.Point(0.690, 1.598), vg.Point(0.990, 3.102), vg.Point(1.148, 3.259)],
                [vg.Point(-2.738, 0.775), vg.Point(-2.738, 1.458), vg.Point(-2.660, 1.850), vg.Point(-2.509, 2.303), vg.Point(-2.271, 2.778), vg.Point(-1.719, 3.331), vg.Point(-1.719, 2.540), vg.Point(-1.790, 1.688), vg.Point(-1.861, 0.903), vg.Point(-2.046, -0.022), vg.Point(-2.525, -0.022), vg.Point(-2.660, 0.383), vg.Point(-2.738, 0.775)]]
    start = time.perf_counter()
    g = vg.VisGraph()
    g.build(polys)
    path = g.shortest_path(vg.Point(0.430, 1.068), vg.Point(7.767, -9.629))
    elapsed_time = (time.perf_counter() - start) * 1e6
    print(f"Calc time: {elapsed_time:.2f} µs")
    print(path)

example_vg2()