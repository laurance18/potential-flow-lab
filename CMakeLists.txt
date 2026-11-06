cmake_minimum_required(VERSION 3.25)

project(PotentialFlowLab
    VERSION 0.3.0
    DESCRIPTION "Interactive 2D potential-flow simulator"
    LANGUAGES C CXX)

if(NOT MSVC)
    message(FATAL_ERROR "Potential Flow Lab is set up for MSVC (Visual Studio 2022) only.")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Standalone exe: static CRT for us and every dependency (must be set before FetchContent).
set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build all dependencies as static libraries" FORCE)

# Outputs land in build/<preset>/bin/<Config>/.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

option(PFL_SHADERS_FROM_DISK "Debug builds read shaders from the source tree (hot reload) before falling back to embedded copies" ON)
option(PFL_BUILD_TESTS "Build the doctest unit tests" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CompilerWarnings)
include(EmbedResources)
include(Dependencies)

add_subdirectory(src/core)
add_subdirectory(src/render)
add_subdirectory(src/ui)
add_subdirectory(src/io)
add_subdirectory(src/app)

if(PFL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

include(Dist)
