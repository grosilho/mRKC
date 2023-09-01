cmake_minimum_required(VERSION 3.17...3.24)

file(MAKE_DIRECTORY build)

project(MultirateIntegrators
        VERSION 1.0 
        DESCRIPTION "Implementatiion of explicit stabilized multirate integrators adn other more classical methods." 
        LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_BUILD_TYPE "Release")

option(BUILD_SHARED_LIBS "Build shared (ON) or static (OFF) libraries" OFF)

add_executable(${PROJECT_NAME} src/main.cpp)

#download submodules
find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    else()
    endif()
endif()

find_package(Eigen3 REQUIRED CONFIG) # Eigen3::Eigen (v 3.4.1 is compatible)

add_subdirectory(src)
target_include_directories(${PROJECT_NAME} 
    PUBLIC include
    PUBLIC external/eigen
    PUBLIC external/getpot
)
target_link_directories(${PROJECT_NAME} PRIVATE src)
target_link_libraries(${PROJECT_NAME} PUBLIC MainHeader parameters)


install(TARGETS ${PROJECT_NAME} DESTINATION .)


