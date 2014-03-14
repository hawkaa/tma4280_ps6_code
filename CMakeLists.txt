project(poisson)

cmake_minimum_required(VERSION 2.6)

# We use the C language
enable_language(C)
enable_language(CXX)
enable_language(Fortran)

# CMake Options to enable/disable MPI and OpenMP
option(ENABLE_OPENMP "Enable OpenMP support?" ON)
option(ENABLE_MPI    "Enable MPI support?" ON)

if(ENABLE_MPI)
  find_package(MPI)
endif()

if(MPI_FOUND)
  add_definitions(-DHAVE_MPI=1)
  set(INCLUDES ${INCLUDES} ${MPI_INCLUDE_PATH})
  set(DEPLIBS_C ${DEPLIBS_C} ${MPI_C_LIBRARIES})
  set(DEPLIBS_CXX ${DEPLIBS_CXX} ${MPI_CXX_LIBRARIES})
  set(DEPLIBS_F ${DEPLIBS_F} ${MPI_Fortran_LIBRARIES})
endif()

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OPENMP_FOUND)
    add_definitions(-DHAVE_OPENMP=1)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  endif()
endif()

# Add gtest
add_subdirectory(gtest)

# Send the include directories to the compiler
include_directories(${INCLUDES} include)

# Add libraries
add_library(ps6_common_library ps6_common_library.c)
add_library(fst fst.f)

# The C based poisson application
add_executable(poisson poisson.c)

# The convergence test program
add_executable(convergence convergence.c)

# The matrix transpose test program
add_executable(transposetest transposetest.cpp)


# Link to libraries
target_link_libraries(ps6_common_library fst ${DEPLIBS_C})
target_link_libraries(poisson ps6_common_library ${DEPLIBS_C})
target_link_libraries(convergence ps6_common_library ${DEPLIBS_C})

target_link_libraries(transposetest ps6_common_library ${DEPLIBS_CXX} gtest gtest_main pthread)
