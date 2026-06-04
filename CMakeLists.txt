cmake_minimum_required(VERSION 3.14)

set(target_name HygroThermFEM)

cmake_policy(SET CMP0048 NEW)

project(${target_name} VERSION 1.0.61 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)

# Local development mode. The `local` preset sets LBNL_LOCAL_SIBLINGS=ON (a cache
# variable, so it cascades into FetchContent sub-builds). When it is on, prefer a
# sibling working copy (../<Name>) of each of THIS project's direct dependencies,
# so the whole dependency graph can be built from local disk during development.
# A dependency whose sibling is not checked out falls back to its declared remote,
# independently. Propagation: a sibling dependency that also honors this flag wires
# up its own siblings in turn. When the flag is off (e.g. the `default` preset),
# nothing here runs and every dependency comes from its remote.
# Note: each pair is "<folder>=<FETCHCONTENT var suffix>"; they differ for eigen
# (folder) vs EIGENEXTERNAL (its FetchContent_Declare name, uppercased).
if(LBNL_LOCAL_SIBLINGS)
    foreach(_sibling "KeffCavity=KEFFCAVITY" "googletest=GOOGLETEST" "eigen=EIGENEXTERNAL")
        string(REPLACE "=" ";" _sibling "${_sibling}")
        list(GET _sibling 0 _sibling_dir)
        list(GET _sibling 1 _sibling_var)
        if(NOT FETCHCONTENT_SOURCE_DIR_${_sibling_var}
           AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../${_sibling_dir}/CMakeLists.txt")
            set(FETCHCONTENT_SOURCE_DIR_${_sibling_var}
                "${CMAKE_CURRENT_SOURCE_DIR}/../${_sibling_dir}")
            message(STATUS "[local] ${_sibling_dir} -> ${FETCHCONTENT_SOURCE_DIR_${_sibling_var}}")
        endif()
    endforeach()
endif()

option(BUILD_FOR_CODE_COVERAGE OFF)

if (BUILD_FOR_CODE_COVERAGE)
    if (CMAKE_COMPILER_IS_GNUCXX)
        set(CMAKE_CXX_FLAGS "--coverage")
    endif()
    if (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
        set(CMAKE_CXX_FLAGS "-fprofile-instr-generate -fcoverage-mapping")
    endif()
endif()

option(BUILD_FOR_PROFILING OFF)

if(BUILD_FOR_PROFILING)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -O0")
endif()

include(FetchContent)
include(cmake/HTProjectMacros.cmake)
include(cmake/HTinternal_utils.cmake)
include(cmake/HTCompilerFlags.cmake)

# Include Eigen. Declared from its remote; in `local` mode the foreach above sets
# FETCHCONTENT_SOURCE_DIR_EIGENEXTERNAL to the ../eigen sibling, so it is sourced
# from local disk instead. Use ${eigenexternal_SOURCE_DIR} for the include path so
# it resolves wherever the source came from (remote _deps clone or local sibling).
set(Eigen_Branch "5.0.1")
if(NOT TARGET eigen)
    FetchContent_Declare(
            EigenExternal
            GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
            GIT_TAG ${Eigen_Branch}
    )
    FetchContent_MakeAvailable(EigenExternal)
    include_directories(${eigenexternal_SOURCE_DIR})
endif()

# Include KeffCavity
set(BUILD_KeffC_TESTING OFF CACHE BOOL "Build KeffCavity tests" FORCE)
set(KeffCavity_Branch "Version_1.1.3")

if(NOT TARGET KeffCavity)
    FetchContent_Declare(
            KeffCavity
            GIT_REPOSITORY https://github.com/LBNL-ETA/KeffCavity.git
            GIT_TAG ${KeffCavity_Branch}
    )
    FetchContent_MakeAvailable(KeffCavity)
endif()

get_property(WindowsCalcEngine_Branch TARGET KeffCavity PROPERTY WindowsCalcEngine_Branch)
configure_file(doc/HygroThermFEM.md.in ${CMAKE_CURRENT_SOURCE_DIR}/doc/HygroThermFEM.md @ONLY)

# google test will not be created by this project if this is not main project. Main project is expected to handle google test.
get_directory_property(hasParent PARENT_DIRECTORY)

if (hasParent)
    set(BUILD_HTFEM_TESTING OFF)
else ()
    option(BUILD_HTFEM_TESTING "Build testing targets" ON)
endif ()

# add_subdirectory( src )

file(GLOB SOURCES_CXX "src/*.cxx")
file(GLOB SOURCES_HXX "src/*.hxx")
file(GLOB INCLUDE_HXX "include/*.hxx")
list(APPEND SOURCES ${SOURCES_HXX} ${INCLUDE_HXX} ${SOURCES_CXX})

if (MSVC)
    add_definitions("/wd4996")
    add_compile_definitions(STL_MULTITHREADING)
endif ()

add_library(${target_name} STATIC ${SOURCES})

set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries( ${target_name} PUBLIC KeffCavity)

if (BUILD_HTFEM_TESTING)
    include_directories(include)
    include_directories(tst/helper)
    
    # Prevent GoogleTest from overriding our compiler/linker options
    # when building with Visual Studio
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    if(NOT TARGET gtest)
        FetchContent_Declare(
            googletest
            GIT_REPOSITORY https://github.com/google/googletest.git
            GIT_TAG v1.16.0
        )
        FetchContent_MakeAvailable(googletest)
    endif()    

    # Note: tst/main.cxx is listed literally rather than via GLOB_RECURSE because
    # GLOB_RECURSE interprets "tst/main.cxx" as a path-suffix match and also
    # matches tst/sweep/main.cxx, pulling the sweep entry point into the unit
    # test target and causing LNK2005 (main already defined).
    file(GLOB_RECURSE all_test_src RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "tst/units/*.cxx" "tst/helper/*.hxx" "tst/helper/*.cxx")
    LIST(APPEND test_src "tst/main.cxx" ${all_test_src})

    CREATE_TEST_TARGETS_HTFEM(${target_name} "${test_src}" ${target_name})

    #add_executable(${target_name}_test ${test_src})

    #target_link_libraries(${target_name}_test gtest ${target_name})

    # Solver parameter sweep — separate executable, not part of normal test suite
    file(GLOB_RECURSE sweep_src RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "tst/sweep/*.cxx" "tst/helper/*.hxx" "tst/helper/*.cxx")
    add_executable(${target_name}_sweep ${sweep_src})
    target_include_directories(${target_name}_sweep PRIVATE tst/helper)
    target_link_libraries(${target_name}_sweep ${target_name} gtest)
endif ()

config_compiler_and_linker_HTFEM()

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
endif()