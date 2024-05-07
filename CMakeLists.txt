cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(hello VERSION 0.0.0)

configure_file(helloConfig.h.in helloConfig.h)
#string(TIMESTAMP COMPILE_TIME %Y%m%d-%H%M%S)

# set the project name
project(helloCFD)


#add_subdirectory(libMeshb)
#list(APPEND EXTRA_LIBS libMeshb)


file(GLOB_RECURSE SRCS CONFIGURE_DEPENDS src/*.cpp)
include_directories(
        ${PROJECT_SOURCE_DIR}/include/
)

# add the executable
add_executable(${PROJECT_NAME} ${SRCS})

target_link_libraries(${PROJECT_NAME} PUBLIC cgns)
target_link_libraries ( ${PROJECT_NAME} PUBLIC ${EXTRA_LIBS})
target_include_directories(${PROJECT_NAME} PUBLIC
                           ${PROJECT_BINARY_DIR}
                           )




