cmake_minimum_required(VERSION 3.16)
project(mesher)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_FIND_FRAMEWORK LAST) # under some circumstances, gdal can be in the frameworks dir that super breaks everything
endif()

option(BUILD_CGAL "Build CGAL and patch to reproduce mesher paper"  OFF )

set(CMAKE_BUILD_TYPE "Release")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

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
        src-cxx/version.h.in
        src-cxx/version.h
)

LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}")

set(EXTERNAL_LINK "")

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
    find_package(CGAL REQUIRED )
    list(APPEND EXTERNAL_LINK CGAL::CGAL)
endif()

# https://www.cgal.org/2024/10/22/cgal601/
if(${CGAL_MAJOR_VERSION} LESS 6)
    # https://www.cgal.org/2013/04/11/cgal-42/
    include(${CGAL_USE_FILE})
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


find_package(GDAL 3.0 REQUIRED)
list(APPEND EXTERNAL_LINK GDAL::GDAL)

find_package(Boost
            COMPONENTS
            program_options
            filesystem
            REQUIRED)

list(APPEND EXTERNAL_LINK Boost::boost Boost::filesystem Boost::program_options)

# we will need to be able to tell the python metis where the shared objects live
set(metis_LIBS "")

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")
find_package(metis 5.1 REQUIRED)
set(metis_LIBS "${METIS_LIBRARY}")



file(WRITE "${CMAKE_BINARY_DIR}/build-bin/metis-config.py" "METIS_DLL=\'${metis_LIBS}\'")

set(SOURCE_FILES src-cxx/mesher.cpp src-cxx/triangle.cpp src-cxx/raster.cpp)


add_executable(
         meshercxx
         ${SOURCE_FILES}
)
set_target_properties(
        meshercxx
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build-bin"
)
target_compile_features(
        meshercxx
        PRIVATE
        cxx_std_14
)
target_include_directories(
        meshercxx
        PRIVATE
        ${CMAKE_BINARY_DIR}/src-cxx # for clion generated files / out of source builds
)
target_link_libraries(
        meshercxx
        m
        ${EXTERNAL_LINK}
)

install(TARGETS meshercxx RUNTIME)
install(FILES ${CMAKE_BINARY_DIR}/build-bin/metis-config.py DESTINATION bin)
