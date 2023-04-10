cmake_minimum_required(VERSION 3.19)
project(GeometricProcessing_OpenMesh)

set(CMAKE_CXX_STANDARD 14)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    message(STATUS "Compiler is AppleClang")
    set(OPENMESH_PATH "/opt/homebrew/Cellar/open-mesh/9.0/include") # open-mesh path
    set(OPENMESH_LIB_PATH "/opt/homebrew/Cellar/open-mesh/9.0/lib")
    set(OMP_PATH  "/opt/homebrew/Cellar/libomp/15.0.7/include") # openmp path
    set(EIGEN_PATH  "/Users/zss/3rdParty/eigen") # eigen path
    include_directories(${OPENMESH_PATH})
    link_directories(${OPENMESH_LIB_PATH})
    include_directories(${OMP_PATH})
    include_directories(${EIGEN_PATH})

endif()

file (GLOB Algorithm_File
        src/*
        )

add_executable(GeometricProcessing_OpenMesh main.cpp
        ${Algorithm_File}
        )

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    target_link_libraries(${CMAKE_PROJECT_NAME}
            optimized OpenMeshCore
            optimized OpenMeshTools
            )


endif()
