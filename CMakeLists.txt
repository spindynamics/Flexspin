cmake_minimum_required(VERSION 3.26)

# C is in the language list only so a private json-c (a C library) can be compiled
# when no system copy is installed. Flexspin itself is pure C++.
project(Flexspin LANGUAGES C CXX)

# ---------------------------------------------------------------------------
# Language level
# ---------------------------------------------------------------------------
# The physics sources are C++11. Set the standard through target properties
# rather than CMAKE_CXX_FLAGS so it travels with the target, and keep
# CXX_EXTENSIONS off so we get "-std=c++11" and not "-std=gnu++11".
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ---------------------------------------------------------------------------
# Build type
# ---------------------------------------------------------------------------
# The old Makefile compiled everything with an unconditional -O3. Default to
# Release (which is -O3 for GCC/Clang) so a bare "cmake -S . -B build" is as fast
# as the build it replaces, while still letting the user override it.
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
                 Debug Release RelWithDebInfo MinSizeRel)
endif()

option(FLEXSPIN_WARNINGS "Compile with -Wall -Wextra" OFF)
option(FLEXSPIN_INSTALL  "Generate install rules for Flexspin" ON)

set(FLEXSPIN_JSONC "auto" CACHE STRING "How to obtain json-c: auto, system, or vendored")
set_property(CACHE FLEXSPIN_JSONC PROPERTY STRINGS auto system vendored)
set(FLEXSPIN_JSONC_GIT_TAG "json-c-0.19-20260627"
    CACHE STRING "git tag of json-c to compile when vendoring")
if(NOT FLEXSPIN_JSONC MATCHES "^(auto|system|vendored)$")
    message(FATAL_ERROR
            "FLEXSPIN_JSONC must be one of: auto, system, vendored (got '${FLEXSPIN_JSONC}')")
endif()

# ---------------------------------------------------------------------------
# MPI
# ---------------------------------------------------------------------------
# find_package(MPI) replaces the old Makefile's hardcoded "mpicxx". The
# MPI::MPI_CXX imported target carries the right include directories, compile
# options and libraries whether or not CMAKE_CXX_COMPILER is itself a wrapper.
find_package(MPI REQUIRED COMPONENTS CXX)

# ---------------------------------------------------------------------------
# json-c
# ---------------------------------------------------------------------------
# CMake ships no FindJSONC module, so there are three routes to the same library:
#
#   1. system config package  - json-c >= 0.15 built with its own CMake installs
#      json-c-config.cmake exporting json-c::json-c.
#   2. system headers + lib   - what older distro packages ship (RHEL 8's
#      json-c-devel, json-c 0.13/0.14): headers and json-c.pc, no config file.
#   3. vendored source        - compile json-c from a pinned git tag into the
#      build tree, so a machine with no json-c package at all still builds.
#
# "auto" tries 1, then 2, then 3. The package name in route 1 must be given in
# lower case: config-package lookup is case-sensitive and json-c installs
# "json-c-config.cmake".
set(FLEXSPIN_JSONC_TARGET "")
set(FLEXSPIN_JSONC_INCLUDES "")

if(NOT FLEXSPIN_JSONC STREQUAL "vendored")
    find_package(json-c CONFIG QUIET)
endif()

if(TARGET json-c::json-c)
    set(FLEXSPIN_JSONC_TARGET json-c::json-c)
    message(STATUS "Flexspin: json-c via its CMake package")

elseif(NOT FLEXSPIN_JSONC STREQUAL "vendored")
    find_path(JSONC_INCLUDE_DIR NAMES json-c/json.h)
    find_library(JSONC_LIBRARY NAMES json-c)
    if(JSONC_INCLUDE_DIR AND JSONC_LIBRARY)
        add_library(json-c::json-c INTERFACE IMPORTED)
        set_target_properties(json-c::json-c PROPERTIES
            INTERFACE_INCLUDE_DIRECTORIES "${JSONC_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "${JSONC_LIBRARY}")
        set(FLEXSPIN_JSONC_TARGET json-c::json-c)
        message(STATUS "Flexspin: json-c ${JSONC_LIBRARY} (${JSONC_INCLUDE_DIR})")
    endif()
endif()

if(NOT FLEXSPIN_JSONC_TARGET)
    if(FLEXSPIN_JSONC STREQUAL "system")
        message(FATAL_ERROR
                "FLEXSPIN_JSONC=${FLEXSPIN_JSONC}, but no system json-c was found.\n"
                "The runtime library alone is not enough: every source includes objects.hpp,\n"
                "which includes <json-c/json.h>. Install the development files\n"
                "  Debian/Ubuntu : sudo apt install libjson-c-dev\n"
                "  RHEL/Fedora   : sudo dnf install json-c-devel\n"
                "or drop -DFLEXSPIN_JSONC=auto to compile a private copy from source.")
    endif()

    include(FetchContent)

    # Steer json-c's own build. These are cache variables because json-c declares
    # them with option() and we want our choice to win over the default.
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)  # link it in: no runtime .so lookup
    set(BUILD_STATIC_LIBS OFF CACHE BOOL "" FORCE)   # ...and skip the duplicate
                                                     # json-c-static target
    set(DISABLE_WERROR ON CACHE BOOL "" FORCE)       # a dependency's -Werror must
                                                     # never be able to fail our build
    set(BUILD_APPS OFF CACHE BOOL "" FORCE)          # json-c's own CLI tools
    FetchContent_Declare(json-c
        GIT_REPOSITORY https://github.com/json-c/json-c.git
        GIT_TAG         "${FLEXSPIN_JSONC_GIT_TAG}"
        GIT_SHALLOW     TRUE
    )
    FetchContent_MakeAvailable(json-c)

    # json-c generates json.h into its *build* tree but includes its siblings
    # ("arraylist.h", "json_object.h", ...) by plain quoted name, and those live in
    # the *source* tree. So neither tree alone satisfies the <json-c/json.h>
    # include that objects.hpp uses -- only the installed layout does, with every
    # public header under a json-c/ subdirectory. Assemble exactly that layout in
    # the build tree with symlinks (copies where symlinks are unavailable), so
    # objects.hpp stays identical whether json-c came from the system or here.
    set(_jsonc_stage "${CMAKE_CURRENT_BINARY_DIR}/json-c-vendored/include/json-c")
    file(MAKE_DIRECTORY "${_jsonc_stage}")
    foreach(_dir IN ITEMS "${json-c_SOURCE_DIR}" "${json-c_BINARY_DIR}")
        file(GLOB _headers "${_dir}/*.h")
        foreach(_header IN LISTS _headers)
            get_filename_component(_name "${_header}" NAME)
            file(CREATE_LINK "${_header}" "${_jsonc_stage}/${_name}"
                 SYMBOLIC COPY_ON_ERROR)
        endforeach()
    endforeach()

    set(FLEXSPIN_JSONC_TARGET json-c)
    set(FLEXSPIN_JSONC_INCLUDES "${CMAKE_CURRENT_BINARY_DIR}/json-c-vendored/include")
    message(STATUS "Flexspin: json-c ${FLEXSPIN_JSONC_GIT_TAG} compiled from source "
                   "(${json-c_BINARY_DIR})")
endif()

# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------
# Listed explicitly rather than with file(GLOB): CMake does not re-run configure
# when a GLOB pattern starts matching new files, so a glob silently omits newly
# added .cpp files until somebody happens to touch this file. Add sources here.
add_executable(Flexspin
    src/CalcHeff.cpp
    src/LLG_Heun1.cpp
    src/LLG_Heun2.cpp
    src/LLG_RK4.cpp
    src/LLG_Symplectic.cpp
    src/Ndip.cpp
    src/Rand_dist.cpp
    src/engine.cpp
    src/evol_filename.cpp
    src/junction_init.cpp
    src/main.cpp
    src/path_resolver.cpp
    src/prints.cpp
    src/pulse.cpp
    src/read_input.cpp
    src/simulation_range.cpp
)

# objects.hpp and printing.hpp both live in src/.
target_include_directories(Flexspin PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/src")

# Empty unless json-c was vendored; a system json-c carries its own includes.
if(FLEXSPIN_JSONC_INCLUDES)
    target_include_directories(Flexspin PRIVATE ${FLEXSPIN_JSONC_INCLUDES})
endif()

target_link_libraries(Flexspin PRIVATE MPI::MPI_CXX ${FLEXSPIN_JSONC_TARGET})

# sqrt/sin/floor/... on the platforms that still keep libm separate from libc.
if(UNIX)
    target_link_libraries(Flexspin PRIVATE m)
endif()

# Off by default: the sources still carry unused parameters (e.g. tot_layers in
# the LLG solvers), so -Wextra is noisy until those are cleaned up.
if(FLEXSPIN_WARNINGS)
    target_compile_options(Flexspin PRIVATE -Wall -Wextra)
endif()

if(FLEXSPIN_INSTALL)
    include(GNUInstallDirs)
    install(TARGETS Flexspin RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
endif()