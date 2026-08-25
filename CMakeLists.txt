cmake_minimum_required(VERSION 3.14)

set(target_name HygroThermFEM)

cmake_policy(SET CMP0048 NEW)

project(${target_name} VERSION 1.0.74 LANGUAGES CXX)

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
    foreach(_sibling "KeffCavity=KEFFCAVITY" "LBNLCPPCommon=LBNLCPPCOMMON" "googletest=GOOGLETEST" "eigen=EIGENEXTERNAL" "FileParse=FILEPARSE" "xmlParser=XMLPARSER" "json=NLOHMANN_JSON")
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

option(BUILD_FOR_CODE_COVERAGE "Build with code-coverage instrumentation" OFF)

# Append rather than assign: a plain set() here dropped whatever the preset, the toolchain
# file or CMAKE_BUILD_TYPE had already put in CMAKE_CXX_FLAGS. The linker needs the same
# flag as the compiler, otherwise the instrumented objects fail to link against gcov/profile.
# MATCHES "Clang" covers both Clang and AppleClang; plain Clang was not handled before.
if (BUILD_FOR_CODE_COVERAGE)
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --coverage")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-instr-generate -fcoverage-mapping")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fprofile-instr-generate")
        set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fprofile-instr-generate")
    else()
        message(WARNING "BUILD_FOR_CODE_COVERAGE is ON but ${CMAKE_CXX_COMPILER_ID} has no "
                        "coverage flags configured; building without instrumentation.")
    endif()
endif()

option(BUILD_FOR_PROFILING "Build with profiling-friendly optimization settings" OFF)

if(BUILD_FOR_PROFILING)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
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

# Include KeffCavity. It also owns the Windows-CalcEngine pin, re-exported as a target
# property and read back below, so the version is chosen in exactly one place.
set(BUILD_KeffC_TESTING OFF CACHE BOOL "Build KeffCavity tests" FORCE)
# Branch pin while the Windows/Clang -fPIC fix is under test. KeffCavity carries
# the same fix and re-exports the Windows-CalcEngine pin. Re-pin to a released
# tag once it lands.
set(KeffCavity_Branch "Version_1.1.7")

if(NOT TARGET KeffCavity)
    FetchContent_Declare(
            KeffCavity
            GIT_REPOSITORY https://github.com/LBNL-ETA/KeffCavity.git
            GIT_TAG ${KeffCavity_Branch}
    )
    FetchContent_MakeAvailable(KeffCavity)
endif()

# LBNLCPPCommon: ExpectedExt/OptionalExt etc. (the C++20 stand-ins for std::expected/
# std::optional) used by the errorest component. Header-only; no C++ standard bump.
set(LBNLCPPCommon_Branch "v0.19")
if(NOT TARGET LBNLCPPCommon)
    FetchContent_Declare(
            LBNLCPPCommon
            GIT_REPOSITORY https://github.com/LBNL-ETA/LBNLCPPCommon.git
            GIT_TAG ${LBNLCPPCommon_Branch}
    )
    FetchContent_MakeAvailable(LBNLCPPCommon)
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

file(GLOB_RECURSE SOURCES_CXX CONFIGURE_DEPENDS "src/*.cxx")
file(GLOB_RECURSE SOURCES_HXX CONFIGURE_DEPENDS "src/*.hxx")
file(GLOB INCLUDE_HXX CONFIGURE_DEPENDS "include/*.hxx")
list(APPEND SOURCES ${SOURCES_HXX} ${INCLUDE_HXX} ${SOURCES_CXX})

if (MSVC)
    add_definitions("/wd4996")
    add_compile_definitions(STL_MULTITHREADING)
endif ()

add_library(${target_name} STATIC ${SOURCES})

set_target_properties(${target_name} PROPERTIES LINKER_LANGUAGE CXX)

target_link_libraries( ${target_name} PUBLIC KeffCavity LBNLCPPCommon)

# src/ is grouped into component subfolders; put each on the include path so the
# existing flat #includes ("Foo.hxx") keep resolving across components. (Step 2 of
# the structure refactor; the curated public facades come in Step 3.)
set(HTFEM_COMPONENT_DIRS
    ${CMAKE_CURRENT_SOURCE_DIR}/src/boundary
    ${CMAKE_CURRENT_SOURCE_DIR}/src/common
    ${CMAKE_CURRENT_SOURCE_DIR}/src/domain
    ${CMAKE_CURRENT_SOURCE_DIR}/src/elements
    ${CMAKE_CURRENT_SOURCE_DIR}/src/materials
    ${CMAKE_CURRENT_SOURCE_DIR}/src/math
    ${CMAKE_CURRENT_SOURCE_DIR}/src/mesh
    ${CMAKE_CURRENT_SOURCE_DIR}/src/simulation
    ${CMAKE_CURRENT_SOURCE_DIR}/src/errorest)
# Public header surface carried by the target itself, so consumers get it by linking:
#   include/          the HygroThermFEM2D.hxx umbrella and the curated HygroThermFEM/ facades
#   src/<component>/  the flat historical includes those facades resolve to
#   Eigen             SquareMatrix.hxx includes <Eigen/SparseCore>, so Eigen is part of this
#                     target's public headers and has to travel with them
# Previously include/ was added with directory scope only, and only under BUILD_HTFEM_TESTING,
# while Eigen was added with directory scope at all times. Neither reached a consumer, so every
# consumer re-derived both paths by hand (Windows-Tools does exactly that in Products/Therm and
# Products/Window). Those hand-wired entries stay harmless once this exports the same paths.
FetchContent_GetProperties(EigenExternal)

target_include_directories(${target_name} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${HTFEM_COMPONENT_DIRS})

# SYSTEM so Eigen's own warnings do not land in consumers' builds.
if(eigenexternal_SOURCE_DIR)
    target_include_directories(${target_name} SYSTEM PUBLIC ${eigenexternal_SOURCE_DIR})
endif()

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
            GIT_TAG v1.17.0
        )
        FetchContent_MakeAvailable(googletest)
    endif()

    # FileParse: the house serialization library (JSON/XML behind one declarative interface).
    # Test-only dependency, used solely to load the error-estimator regression fixtures
    # (tst/units/errorest/fixtures/*.json) via Common::loadFromFile. It never reaches the
    # HygroThermFEM library, so consumers (BUILD_HTFEM_TESTING off) never pull it.
    # FileParse PUBLIC-links xmlParser + nlohmann_json; in local mode the foreach above
    # points all three at their ../ siblings.
    set(FileParse_Branch "Version_1.1.6")
    set(BUILD_FileParse_tests OFF CACHE BOOL "Build FileParse tests" FORCE)
    if(NOT TARGET FileParse)
        FetchContent_Declare(
            FileParse
            GIT_REPOSITORY https://github.com/LBNL-ETA/FileParse.git
            GIT_TAG ${FileParse_Branch}
        )
        FetchContent_MakeAvailable(FileParse)
    endif()

    # Note: tst/main.cxx is listed literally rather than via GLOB_RECURSE because
    # GLOB_RECURSE interprets "tst/main.cxx" as a path-suffix match and also
    # matches tst/sweep/main.cxx, pulling the sweep entry point into the unit
    # test target and causing LNK2005 (main already defined).
    file(GLOB_RECURSE all_test_src RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "tst/units/*.cxx" "tst/helper/*.hxx" "tst/helper/*.cxx")
    LIST(APPEND test_src "tst/main.cxx" ${all_test_src})

    CREATE_TEST_TARGETS_HTFEM(${target_name} "${test_src}" ${target_name})

    # errorest regression fixtures are loaded at runtime via FileParse (test-only).
    target_link_libraries(${target_name}_tests FileParse)
    target_compile_definitions(${target_name}_tests PRIVATE
        ERROREST_FIXTURE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tst/units/errorest/fixtures")

    # Register the suite with CTest so `ctest --preset <name>` and IDE test runners see it.
    # ADD_GOOGLE_TESTS_HTFEM above cannot do this: it scrapes TEST()/TEST_F() out of files
    # matching \.(cpp|cc)$, and every test here is .cxx, so it always registered nothing --
    # and without enable_testing() its add_test() calls would have been dropped anyway.
    # One entry running the whole gtest binary, matching the HygroTherm sibling's setup.
    enable_testing()
    add_test(NAME ${target_name}_tests COMMAND ${target_name}_tests)

    #add_executable(${target_name}_test ${test_src})

    #target_link_libraries(${target_name}_test gtest ${target_name})

    # Solver parameter sweep — separate executable, not part of normal test suite
    file(GLOB_RECURSE sweep_src RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "tst/sweep/*.cxx" "tst/helper/*.hxx" "tst/helper/*.cxx")
    add_executable(${target_name}_sweep ${sweep_src})
    target_include_directories(${target_name}_sweep PRIVATE tst/helper)
    target_link_libraries(${target_name}_sweep ${target_name} gtest)

    # --- Header self-containment check -----------------------------------
    # Every src/ header must compile on its own (i.e. include everything it
    # uses). One generated TU per header just #includes it; building this
    # object library fails if any header is missing an include. It links the
    # library only to inherit its include directories (Eigen, KeffCavity, ...).
    set(_sc_dir ${CMAKE_BINARY_DIR}/selfcontained)
    file(MAKE_DIRECTORY ${_sc_dir})
    file(GLOB_RECURSE PUBLIC_HXX "include/*.hxx")
    set(_sc_tus "")
    foreach(_hdr ${SOURCES_HXX} ${PUBLIC_HXX})
        get_filename_component(_hname ${_hdr} NAME_WE)
        set(_tu ${_sc_dir}/sc_${_hname}.cxx)
        file(WRITE ${_tu} "#include \"${_hdr}\"\n")
        list(APPEND _sc_tus ${_tu})
    endforeach()
    add_library(${target_name}-selfcontained OBJECT ${_sc_tus})
    target_compile_features(${target_name}-selfcontained PRIVATE cxx_std_20)
    # Include dirs only (no link) so this compiles headers without building any lib.
    target_include_directories(${target_name}-selfcontained PRIVATE
        ${HTFEM_COMPONENT_DIRS}
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/include/HygroThermFEM
        ${eigenexternal_SOURCE_DIR}
        $<TARGET_PROPERTY:KeffCavity,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:LBNLCPPCommon,INTERFACE_INCLUDE_DIRECTORIES>
        $<TARGET_PROPERTY:Windows-CalcEngine,INTERFACE_INCLUDE_DIRECTORIES>)
endif ()

# Link-time optimization. The innermost element-assembly loop calls SquareMatrix accessors
# and operators that live in a different translation unit, so those calls cannot be inlined
# without it. Applied to this project's own targets only; the FetchContent dependencies keep
# whatever settings they declare for themselves.
option(HTFEM_ENABLE_LTO "Enable interprocedural optimization for Release builds" ON)

if (HTFEM_ENABLE_LTO)
    include(CheckIPOSupported)
    # LANGUAGES CXX because the project enables only C++; the default check would
    # probe the (unconfigured) C compiler and always report IPO as unavailable.
    check_ipo_supported(RESULT HTFEM_IPO_SUPPORTED OUTPUT HTFEM_IPO_ERROR LANGUAGES CXX)
    if (HTFEM_IPO_SUPPORTED)
        foreach(ipo_target ${target_name} ${target_name}_tests ${target_name}_sweep)
            if (TARGET ${ipo_target})
                set_property(TARGET ${ipo_target}
                             PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
            endif()
        endforeach()
    else()
        message(STATUS "[HTFEM] IPO unavailable, building without LTO: ${HTFEM_IPO_ERROR}")
    endif()
endif()

config_compiler_and_linker_HTFEM()

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")
endif()