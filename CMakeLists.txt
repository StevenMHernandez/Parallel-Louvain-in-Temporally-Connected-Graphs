cmake_minimum_required(VERSION 3.2)
project(LouvainMPI)

# find header & source
file(GLOB_RECURSE SOURCE "main.cpp")

add_executable(${PROJECT_NAME}
        main.cpp
        ${HEADER}
        )


find_package(MPI REQUIRED)

set(CMAKE_C_COMPILER /usr/local/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/local/bin/mpic++)

include_directories(${MPI_INCLUDE_PATH})
target_link_libraries(${PROJECT_NAME} ${MPI_LIBRARIES})

source_group("Header include" FILES ${HEADER})
source_group("Source src"     FILES ${SOURCE})
