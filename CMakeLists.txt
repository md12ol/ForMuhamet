cmake_minimum_required(VERSION 3.30.5)
project(ForMuhamet)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

## Add source files
#add_executable(SDA_test
#        SDA_test.cpp
#        SDA/SDA.cpp
#        constants.h
#)
#
## Add source files
#add_executable(Graph_test
#        Graph_test.cpp
#        Graph/Graph.cpp
#        Combo_test.cpp
#        constants.h
#)

find_package(Pagmo CONFIG REQUIRED)

add_executable(pagmo_test pagmoTest.cpp)
target_link_libraries(pagmo_test PRIVATE Pagmo::pagmo)

add_executable(pagmo_sga
        sga.cpp
        SDA/SDA.cpp
        Graph/Graph.cpp
        rng.h
)
target_link_libraries(pagmo_sga PRIVATE Pagmo::pagmo)

# Das ist dein neues NSGA-II Programm
add_executable(pagmo_nsga2
        nsga2.cpp
        SDA/SDA.cpp
        Graph/Graph.cpp
        rng.h
)
# Pagmo auch hier verlinken!
target_link_libraries(pagmo_nsga2 PRIVATE Pagmo::pagmo)

add_executable(Test
        main.cpp
        Graph/Graph.cpp
        SDA/SDA.cpp
        rng.h
)

