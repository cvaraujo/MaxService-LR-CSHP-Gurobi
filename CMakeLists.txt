cmake_minimum_required(VERSION 3.5)
project(MaxServiceRL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/gurobi900/linux64/include /opt/boost_1_72_0)

link_directories(/opt/gurobi900/linux64/lib /opt/boost_1_72_0)

set(SOURCE_FILES main.cpp src/Graph.cpp headers/Graph.h headers/Include.h src/Model.cpp headers/Model.h src/Arc.cpp headers/Arc.h)

add_executable(MaxServiceRL ${SOURCE_FILES})

target_link_libraries(MaxServiceRL gurobi_c++ gurobi90 pthread m)
