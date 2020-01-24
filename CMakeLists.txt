cmake_minimum_required(VERSION 3.0)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

project(FEM)

add_executable(
    program
    Aggregation.cpp
    Aggregation.h
    Element.cpp
    Element.h
    Grid.cpp
    Grid.h
    GlobalData.cpp
    GlobalData.h
    Jakobian.cpp
    Jakobian.h
    MES.cpp
    Node.cpp
    Node.h
    pch.cpp
    pch.h
    universalElement.cpp
    universalElement.h
)