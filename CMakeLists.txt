cmake_minimum_required(VERSION 2.8)
project(mesher)

option(BUILD_CGAL "Build CGAL and patch to reproduce mesher paper"  OFF )

# lovely CMake script to integrate git hashes
# http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/
# Get the current working branch
# Generate gitrevision.hh if Git is available
# and the .git directory is present
# this is the case when the software is checked out from a Git repo
find_program(GIT_SCM git DOC "Git version control")
mark_as_advanced(GIT_SCM)
find_file(GITDIR NAMES .git PATHS ${CMAKE_SOURCE_DIR} NO_DEFAULT_PATH)
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${GITDIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${GITDIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

configure_file(
        src/version.h.in
        src/version.h
)


#as per http://cgal-discuss.949826.n4.nabble.com/CMake-and-flags-td949906.html
#don't override internal settings
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE )


if(BUILD_CGAL)
    configure_file(${CMAKE_SOURCE_DIR}/cgal.4.11.patch ${CMAKE_BINARY_DIR}/cgal.4.11.patch COPYONLY)

    configure_file(CMakeLists_external.txt.in
            lib/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib )

    set(CGAL_DIR "${CMAKE_BINARY_DIR}/lib/CGAL/lib64/CGAL")
    find_package(CGAL REQUIRED  HINTS ${CGAL_DIR})
else()
    set(CGAL_DIR "")
endif()

find_package(CGAL REQUIRED  HINTS ${CGAL_DIR})

if(CGAL_FOUND)
    message(STATUS "Found CGAL ")
    message(STATUS "${CGAL_INCLUDE_DIRS}")
    include(${CGAL_USE_FILE}) #as per https://www.cgal.org/releases.html release 4.2
endif()

#ignore these two under Clion as CGAL will complain
if(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo OR
        CMAKE_BUILD_TYPE MATCHES MinSizeRel OR
        NOT CMAKE_BUILD_TYPE)

    set(CMAKE_BUILD_TYPE "Release")
endif()

#reset these back
if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 ")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#-fsanitize=address
#add_definitions(-DCGAL_DISABLE_ROUNDING_MATH_CHECK=ON)

find_package(GDAL REQUIRED)
find_package(Boost
        COMPONENTS
        program_options
        filesystem
        REQUIRED)

include_directories(
        ${CMAKE_BINARY_DIR}/src # for clion generated files / out of source builds
        ${GDAL_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        ${CGAL_INCLUDE_DIRS}
        ${CGAL_3RD_PARTY_INCLUDE_DIRS}
        )


set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

set(SOURCE_FILES src/mesher.cpp src/triangle.cpp src/raster.cpp)
add_executable(mesher ${SOURCE_FILES})
target_link_libraries(
        mesher
        m
        ${GDAL_LIBRARY}
        ${Boost_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
)
