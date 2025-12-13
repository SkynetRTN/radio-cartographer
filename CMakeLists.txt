cmake_minimum_required(VERSION 3.0)
project(radio_cartographer)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

file(GLOB SOURCES "src/*.cpp") 

include_directories(include)

add_executable(radio-cartographer ${SOURCES})

target_link_libraries(radio-cartographer fftw3 cfitsio CCfits pthread )

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
FetchContent_MakeAvailable(json)
target_link_libraries(radio-cartographer nlohmann_json::nlohmann_json)
