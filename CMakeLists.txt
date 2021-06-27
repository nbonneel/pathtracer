PROJECT(RayTrace)
cmake_minimum_required(VERSION 3.1)


set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11
MESSAGE(STATUS "C++11 enabled. ")



#---- wxWidgets
find_package(wxWidgets REQUIRED)
include(${wxWidgets_USE_FILE})

#---- OPENMP
if(APPLE)
    message("-- Setting OpenMP flags on MacOs. Assuming `brew install libomp`")
    if(CMAKE_C_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_C_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    endif()

    if(CMAKE_CXX_COMPILER_ID MATCHES "Clang\$")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I/usr/local/opt/libomp/include" CACHE INTERNAL "OpenMP flags for Xcode toolchain.")
        set(OpenMP_CXX_LIB_NAMES "omp" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
        set(OpenMP_omp_LIBRARY "/usr/local/opt/libomp/lib/libomp.dylib" CACHE INTERNAL "OpenMP lib name for Xcode toolchain.")
    endif()
endif()
find_package(OpenMP REQUIRED)


SET(SRCS
  Geometry.cpp
  MERLBRDFRead.cpp
  PointSet.cpp
  Raytracer.cpp
  TriangleMesh.cpp
  Vector.cpp
  hdrloader.cpp
  hdrwriter.cpp
  mainApp.cpp
  utils.cpp
  fluid.cpp
  )

add_executable(rayTracer ${SRCS})
target_link_libraries(rayTracer PRIVATE ${wxWidgets_LIBRARIES}  OpenMP::OpenMP_CXX)
