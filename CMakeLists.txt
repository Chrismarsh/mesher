cmake_minimum_required(VERSION 3.16)
project(mesher)


option(BUILD_CGAL "Build CGAL and patch to reproduce mesher paper"  OFF )
option(USE_CONAN "Build without using conan dependencies" FALSE)

if(USE_CONAN AND NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")

    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/develop/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake")
    include(${CMAKE_BINARY_DIR}/conan.cmake)
endif()


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CONAN_LIB_STD "libstdc++11")
    message(STATUS "Detected gnu compiler, setting compiler.libcxx=${CONAN_LIB_STD}")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    set(CONAN_LIB_STD "libc++")
    message(STATUS "Detected clang compiler, setting compiler.libcxx=${CONAN_LIB_STD}")
endif ()

if(USE_CONAN)
    conan_add_remote(NAME bincrafters
                URL https://api.bintray.com/conan/bincrafters/public-conan
            )
    conan_add_remote(NAME CHM
            URL https://api.bintray.com/conan/chrismarsh/CHM
            )

    conan_cmake_run(CONANFILE conanfile.py
                    SETTINGS compiler.cppstd=14
                    SETTINGS compiler.libcxx=${CONAN_LIB_STD}
                    BASIC_SETUP
                    CMAKE_TARGETS
                    KEEP_RPATHS
                    NO_OUTPUT_DIRS
                    BUILD missing)
endif()

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
    if(NOT USE_CONAN)
        find_package(CGAL REQUIRED )
        list(APPEND EXTERNAL_LINK CGAL::CGAL)
    else()
        find_package(cgal REQUIRED ) #conan convention is lowercase
        list(APPEND EXTERNAL_LINK cgal::cgal)
    endif()
endif()

if(CGAL_FOUND AND NOT USE_CONAN)
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

if(NOT USE_CONAN)
    find_package(GDAL REQUIRED)
    list(APPEND EXTERNAL_LINK GDAL::GDAL)
else()
    find_package(gdal REQUIRED)
    list(APPEND EXTERNAL_LINK gdal::gdal)

endif()

find_package(Boost
            COMPONENTS
            program_options
            filesystem
            REQUIRED)
if(NOT USE_CONAN)
    list(APPEND EXTERNAL_LINK Boost::boost Boost::filesystem Boost::program_options)
else()
    list(APPEND EXTERNAL_LINK Boost::Boost)
endif()


set(SOURCE_FILES src/mesher.cpp src/triangle.cpp src/raster.cpp)

# only need fancy rpaths if we are shipping our libs alongside the binary via conan
if(USE_CONAN)
    # when we run make install all the RPATHs will be stripped out and set to these
    # n.b. $ORIGIN should have its dollar sign escaped with a backslash to have it end up with proper
    # syntax in the final executable.
    # https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
    # !!!  these must be done before targets are defined in add_executable !!!!
    if (APPLE)
        set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
    else()
        set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib")
    endif()
endif()

add_executable(mesher ${SOURCE_FILES})
set_target_properties(mesher
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/build-bin"
        )
target_compile_features(mesher PRIVATE cxx_std_14)
target_include_directories(
        mesher
        PRIVATE
        ${CMAKE_BINARY_DIR}/src # for clion generated files / out of source builds
)
target_link_libraries(
        mesher
        m
        ${EXTERNAL_LINK}
)



install(TARGETS mesher RUNTIME)

# only ship our libs if built with conan
if(USE_CONAN)
    install(DIRECTORY ${CMAKE_BINARY_DIR}/lib/
            DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endif()