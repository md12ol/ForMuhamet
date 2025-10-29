cmake_minimum_required(VERSION 4.0)
project(SDA)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Add source files
add_executable(SDA_test
        SDA_test.cpp
        SDA/SDA.cpp
)

# Add source files
add_executable(Graph_test
        Graph_test.cpp
        Graph/Graph.cpp
        Combo_test.cpp
)
