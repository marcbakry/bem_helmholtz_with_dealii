cmake_minimum_required(VERSION 3.21)
if(WIN32)
    cmake_policy(SET CMP0057 NEW)
endif()

project(1_acoustique_BEM_Dirichlet LANGUAGES CXX)

set(CMAKE_VERBOSE_MAKEFILE FALSE)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

if(WIN32)
    set(DEAL_II_DIR "D:/Documents/Programmation/CPP/libs/libDEALII/lib/cmake/deal.II")
    set(TBB_DIR "C:/Program Files\ (x86)/Intel/oneAPI/tbb/latest/lib/cmake/tbb")
    find_package(deal.II 9.6.0 HINTS ${DEAL_II_DIR})
    find_package(TBB HINTS ${TBB_DIR})

else()
    find_package(deal.II 9.6.0)
    find_package(TBB)
endif()


if(NOT ${deal.II_FOUND})
    message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()
if(NOT ${TBB_FOUND})
    message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif()

deal_ii_initialize_cached_variables()

add_executable(exe_1_acoustique_BEM_Dirichlet ./src/main.cpp)
target_sources(exe_1_acoustique_BEM_Dirichlet PUBLIC ./src/laplace_kernel.cpp ./src/helmholtz_kernel.cpp ./src/acoustic_bem.cpp)
target_include_directories(exe_1_acoustique_BEM_Dirichlet PUBLIC ./include D:/Documents/Programmation/CPP/libs/libopenblas/include)

deal_ii_setup_target(exe_1_acoustique_BEM_Dirichlet)

TARGET_LINK_LIBRARIES(exe_1_acoustique_BEM_Dirichlet TBB::tbb)

set_property(TARGET exe_1_acoustique_BEM_Dirichlet PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded")

if(WIN32)
    add_custom_command(TARGET exe_1_acoustique_BEM_Dirichlet POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy D:/Documents/Programmation/CPP/libs/libopenblas/bin/libopenblas.dll ${1_acoustique_BEM_Dirichlet_BINARY_DIR})
endif()