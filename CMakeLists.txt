cmake_minimum_required(VERSION 3.5)
project(radio_cartographer)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Debug)

option(BUILD_PYTHON_BINDINGS "Build shared library for Python bindings" OFF)
option(ENABLE_PROFILING "Enable profiling with gprof" OFF)
option(ENABLE_NATIVE_OPTIMIZATION "Enable architecture-specific optimizations" OFF)

if(ENABLE_NATIVE_OPTIMIZATION)
    if(MSVC)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /arch:AVX2") # Valid for most modern CPUs, auto-detection on MSVC is harder
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native")
    endif()
endif()

if(ENABLE_PROFILING)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

set(radio_cartographer_sources
    src/Analysis.cpp
    src/BackgroundCUDA.cpp
    src/Cartographer.cpp
    src/Composite.cpp
    src/Debugger.cpp
    src/FourtyParser.cpp
    src/FunctionalForm.cpp
    src/GBParser.cpp
    src/Map.cpp
    src/NonParametric.cpp
    src/OutputFile.cpp
    src/Pixel.cpp
    src/PreProcessor.cpp
    src/Processor.cpp
    src/ProcessorRFI.cpp
    src/ProcessorTS.cpp
    src/ProcessorThetaGap.cpp
    src/RCR.cpp
    src/Scan.cpp
    src/Survey.cpp
    src/Tools.cpp
    src/timer.cpp
)

set(radio_cartographer_library_type STATIC)
if(BUILD_PYTHON_BINDINGS)
    set(radio_cartographer_library_type SHARED)
endif()

add_library(radio_cartographer ${radio_cartographer_library_type} ${radio_cartographer_sources})

if(BUILD_PYTHON_BINDINGS)
    set_target_properties(radio_cartographer PROPERTIES
        POSITION_INDEPENDENT_CODE ON
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON
    )
endif()

target_include_directories(radio_cartographer
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

find_package(nlohmann_json 3.11 QUIET)

if(NOT nlohmann_json_FOUND)
    include(FetchContent)

    FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz)
    FetchContent_MakeAvailable(json)
    set(radio_cartographer_bundled_json TRUE)
endif()

find_package(Threads REQUIRED)

target_link_libraries(radio_cartographer
    PUBLIC
        fftw3
        cfitsio
        CCfits
        Threads::Threads
        nlohmann_json::nlohmann_json
)

add_executable(radio-cartographer src/Source.cpp)
target_link_libraries(radio-cartographer PRIVATE radio_cartographer)

install(TARGETS radio_cartographer radio-cartographer
    EXPORT radio_cartographerTargets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

if(radio_cartographer_bundled_json)
    install(TARGETS nlohmann_json
        EXPORT radio_cartographerTargets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
endif()

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT radio_cartographerTargets
    NAMESPACE radio_cartographer::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/radio_cartographer
)

configure_package_config_file(
    cmake/radio_cartographer-config.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/radio_cartographer-config.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/radio_cartographer
)

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/radio_cartographer-config.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/radio_cartographer
)
