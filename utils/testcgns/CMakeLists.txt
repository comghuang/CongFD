cmake_minimum_required(VERSION 3.0.0)
project(cgnstest VERSION 0.1.0 LANGUAGES C)

include(CTest)
enable_testing()

add_executable(cgnstest write_grid_str.c)
target_link_libraries(cgnstest cgns)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
