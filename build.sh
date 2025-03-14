#!/bin/bash

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

# Get the directory where the script is located
script_dir="$(cd "$(dirname "${BASH_SOURCE}")" && pwd)"

# Define Build Directory
build_dir="$script_dir/build"

# Try to Create the Build Directory
mkdir -p "$build_dir"

# Clean the Build directory
find "$build_dir" -mindepth 1 -delete

# Navigate to the Build directory
cd "$build_dir"

# Run cmake and make
cmake ..
make

# Return to the previous directory
cd -