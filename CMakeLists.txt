cmake_minimum_required(VERSION 3.0)
project(radio_cartographer)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_COMPILER g++)
file(GLOB SOURCES "*.cpp") 

include_directories(.)

add_executable(main ${SOURCES})

target_link_libraries(main fftw3 cfitsio CCfits pthread)
