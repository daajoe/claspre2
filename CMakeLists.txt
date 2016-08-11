cmake_minimum_required(VERSION 2.6)

project(claspre)
add_definitions( -DWITH_THREADS=0 )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(libclasp libprogram_opts lib_exst)
add_subdirectory(libclasp)
add_subdirectory(libprogram_opts)
add_subdirectory(app)
add_subdirectory(lib_exst)


