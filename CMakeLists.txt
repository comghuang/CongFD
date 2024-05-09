cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(congFD VERSION 0.0.0)

configure_file(helloConfig.h.in helloConfig.h)
#string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)


#add_subdirectory(libMeshb)
#list(APPEND EXTRA_LIBS libMeshb)




#add the executable
add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC libcongFD)

add_subdirectory(src)
add_subdirectory(utils)

add_executable(cgnstest cgnsReadTest.cpp)
target_link_libraries(cgnstest PUBLIC libcongFD)


