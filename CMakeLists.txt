cmake_minimum_required(VERSION 3.0)
project(radio_cartographer)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_CXX_COMPILER /usr/bin/c++)
file(GLOB SOURCES "*.cpp") 

include_directories(.)

add_executable(main ${SOURCES})

target_link_libraries(main fftw3 cfitsio CCfits pthread )

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)
target_link_libraries(main nlohmann_json::nlohmann_json)
