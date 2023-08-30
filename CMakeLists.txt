cmake_minimum_required(VERSION 3.17...3.24)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/")
if(APPLE)
	message(STATUS "Platform: MacOSX")
	set(Eigen3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/eigen-3.4.0-MacOSX")
	set(Gmsh_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/gmsh-4.10.5-MacOSX")
	set(FFTW3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/fftw-3.3.10-MacOSX")
else()
	message(STATUS "Platform: Linux")
#	set(Eigen3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/eigen-3.4.0")
#	set(Gmsh_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/gmsh-4.10.5")
#	set(FFTW3_ROOT "${CMAKE_SOURCE_DIR}/../../Libraries/fftw-3.3.10")
#        set(FFTW3_DIR "${CMAKE_SOURCE_DIR}/cmake")
endif()

project(MultirateIntegrators
        VERSION 1.0 
        DESCRIPTION "Implementatiion of explicit stabilized multirate integrators adn other more classical methods." 
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_OPENMP "Enables OpenMP" ON)
option(BUILD_SHARED_LIBS "Build shared (ON) or static (OFF) libraries" OFF)

find_package(Eigen3 REQUIRED CONFIG) # Eigen3::Eigen
find_package(Gmsh REQUIRED) # ${GMSH_LIB}, ${GMSH_INC}
find_package(FFTW3 REQUIRED) #${FFTW3_INCLUDE_DIRS}, ${FFTW3_LIBRARY_DIRS}, ${FFTW3_LIBRARIES}

if(ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        link_libraries(OpenMP::OpenMP_CXX)
    endif()
endif()

include_directories(BEFORE include)

add_subdirectory(src)

add_executable(main src/main.cpp)
target_link_libraries(main PUBLIC MainHeader parameters)

