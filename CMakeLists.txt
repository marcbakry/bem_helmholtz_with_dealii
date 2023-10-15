cmake_minimum_required(VERSION 3.2)
cmake_policy(SET CMP0076 NEW)

project(1_acoustique_BEM_Dirichlet LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

find_package(deal.II 9.6.0)
find_package(TBB)

deal_ii_initialize_cached_variables()

add_executable(exe_1_acoustique_BEM_Dirichlet ./src/main.cpp)
target_sources(exe_1_acoustique_BEM_Dirichlet PUBLIC ./src/laplace_kernel.cpp ./src/helmholtz_kernel.cpp ./src/acoustic_bem.cpp)
target_include_directories(exe_1_acoustique_BEM_Dirichlet PUBLIC ./include)

deal_ii_setup_target(exe_1_acoustique_BEM_Dirichlet)

TARGET_LINK_LIBRARIES(exe_1_acoustique_BEM_Dirichlet TBB::tbb)