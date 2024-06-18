cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(congFD VERSION 0.0.0)

configure_file(helloConfig.h.in helloConfig.h)
set(CMAKE_CXX_STANDARD 23)
#string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)


#add_subdirectory(libMeshb)
#list(APPEND EXTRA_LIBS libMeshb)




#add the executable
add_executable(${PROJECT_NAME} zoneMain.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC libcongFD)

add_subdirectory(src)
add_subdirectory(utils)



