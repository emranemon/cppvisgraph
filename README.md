<!--
The MIT License (MIT)

Copyright (c) 2023 MIEA MD EMON 
https://github.com/emranemon

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
-->

# cppvisgraph - Visibility Graph and Shortest path
Cppvisgraph is a MIT-licensed Cpp library for building visibility graphs from a list of simple obstacle polygons. The visibility graph algorithm (D.T. Lee) runs in O(n^2 log n) time. The shortest path is found using Djikstra's algorithm. This is a C++ translation of the library 
[pyvisgraph](https://github.com/TaipanRex/pyvisgraph/tree/master).
### Environment
Though this is a header only library, And dependencies for building is nothing but C++ Standard Library. 
However, to build python wrapper, pybind11 is required. Check if it already exists or not with the following command.

`python3 -c "import pybind11; print(pybind11.__version__)"`.
#### Requirements
1. C++ 11
2. CMAKE Ver. 3.10.0
3. Python 3
4. pybind11 Ver. 2.10.4
### Use cases
1. Can directly be used in any cpp projects. 
2. Can also be used in python projects with the advantage of faster execution time.
### Usages
Follow the examples written in `vgexamples.cpp` or `vgexamples.py` files under `examples` directory.
#### For C++ projects
Copy the `src/cppvisgraph` dir to your project directory and use accordingly.
#### For python projects
Copy the `build` dir to your project directory and use accordingly.
### Build
Make sure pybind11 is installed. If not, install with the following command.
```
pip install pybind11
```
Navigate to the project dir and execute.
```
./build.sh
```
Grant permission if necessary.
```
chmod -R 777 .
```
### Run Examples
If the Build was successfull, the examples should run.

**c++**
```
./build/vgexample
```
**python**
```
python3 examples/vgexamples.py
```
### References
1. [pyvisgraph](https://github.com/TaipanRex/pyvisgraph/tree/master)