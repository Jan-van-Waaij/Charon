cmake_minimum_required(VERSION 3.14)

# Honour CMAKE_MSVC_RUNTIME_LIBRARY / MSVC_RUNTIME_LIBRARY (static MSVC runtime).
if(POLICY CMP0091)
    cmake_policy(SET CMP0091 NEW)
endif()

project(charon
    VERSION 1.0.0
    DESCRIPTION "CHARON: number of individuals and drift times from environmental DNA (C++ port of Charon.jl)"
    LANGUAGES CXX)

set(CHARON_JULIA_REFERENCE_VERSION "0.3.2")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# A fully static executable is the default on Linux and Windows (MinGW).
# macOS does not support fully static executables; build there with
# -DCHARON_FETCH_DEPS=ON to link zlib statically, so that the binary only
# depends on libraries that ship with macOS.
if(APPLE)
    set(_charon_static_default OFF)
else()
    set(_charon_static_default ON)
endif()
option(CHARON_STATIC "Build a fully static charon executable" ${_charon_static_default})
option(CHARON_NATIVE "Optimise for the build machine's CPU (-march=native); the binary is then not portable" OFF)
option(CHARON_FETCH_DEPS "Download Eigen and zlib instead of using system packages" OFF)
option(CHARON_BUILD_TESTS "Build the unit tests" ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(CharonDependencies)

find_package(Threads REQUIRED)

configure_file(src/charon/version.hpp.in "${CMAKE_CURRENT_BINARY_DIR}/generated/charon/version.hpp" @ONLY)

add_library(charon_core STATIC
    src/charon/data.cpp
    src/charon/io.cpp
    src/charon/model.cpp
    src/charon/priors.cpp
    src/charon/sampler.cpp)
target_include_directories(charon_core PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_BINARY_DIR}/generated")
target_link_libraries(charon_core PUBLIC charon::eigen charon::zlib Threads::Threads)
target_compile_definitions(charon_core PUBLIC EIGEN_MPL2_ONLY)
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9)
    # std::filesystem lives in a separate library before GCC 9.
    target_link_libraries(charon_core PUBLIC stdc++fs)
endif()

if(MSVC)
    target_compile_options(charon_core PUBLIC /W4 /permissive-)
else()
    target_compile_options(charon_core PUBLIC -Wall -Wextra -Wpedantic)
    if(CHARON_NATIVE)
        target_compile_options(charon_core PUBLIC -march=native)
    endif()
endif()

add_executable(charon src/main.cpp)
target_link_libraries(charon PRIVATE charon_core)

if(CHARON_STATIC)
    if(MSVC)
        set_property(TARGET charon PROPERTY MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    else()
        target_link_options(charon PRIVATE -static)
        if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
            # With static glibc, std::thread needs every pthread symbol linked in.
            target_link_libraries(charon PRIVATE -Wl,--whole-archive -lpthread -Wl,--no-whole-archive)
        endif()
    endif()
endif()
if(NOT MSVC AND NOT APPLE)
    # Strip symbols from release binaries.
    target_link_options(charon PRIVATE $<$<CONFIG:Release>:-s>)
endif()

install(TARGETS charon RUNTIME DESTINATION bin)
install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../README.md" "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.md"
        DESTINATION share/doc/charon)

if(CHARON_BUILD_TESTS)
    enable_testing()
    add_executable(charon_tests tests/test_charon.cpp)
    target_link_libraries(charon_tests PRIVATE charon_core)
    target_compile_definitions(charon_tests PRIVATE
        CHARON_TEST_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests/data")
    add_test(NAME unit COMMAND charon_tests)
    add_test(NAME cli_help COMMAND charon --help)
    add_test(NAME cli_dice_smoke
        COMMAND charon --chains 2 --steps 300 --seed 1 --threads 2 --quiet
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/small.dice.csv"
                "${CMAKE_CURRENT_BINARY_DIR}/smoke_dice.csv")
    add_test(NAME cli_basecount_smoke
        COMMAND charon --chains 1 --steps 200 --seed 2 --quiet
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/small_basecount.txt.gz"
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/small_freq.txt.gz"
                "${CMAKE_CURRENT_BINARY_DIR}/smoke_basecount.csv.gz")
    add_test(NAME cli_exactposterior_smoke
        COMMAND charon exactposterior --n 1:3 --tauC 0.1,0.5 --tauA 0.2:0.4:3 --eps 0.01 --quiet
                "${CMAKE_CURRENT_SOURCE_DIR}/tests/data/small.dice.csv"
                "${CMAKE_CURRENT_BINARY_DIR}/smoke_grid.csv")
endif()

set(CPACK_PACKAGE_NAME "charon")
set(CPACK_PACKAGE_VENDOR "CHARON authors")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "${PROJECT_DESCRIPTION}")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/../LICENSE.md")
if(NOT CPACK_SYSTEM_NAME)
    # Package names such as charon-1.0.0-linux-x86_64.tar.gz, as on the releases page.
    string(TOLOWER "${CMAKE_SYSTEM_NAME}" _charon_os)
    if(_charon_os STREQUAL "darwin")
        set(_charon_os "macos")
    endif()
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" _charon_arch)
    if(_charon_arch MATCHES "^(amd64|x64)$")
        set(_charon_arch "x86_64")
    elseif(_charon_arch STREQUAL "arm64" AND NOT APPLE)
        set(_charon_arch "aarch64")
    endif()
    if(APPLE AND CMAKE_OSX_ARCHITECTURES)
        set(_charon_arch "${CMAKE_OSX_ARCHITECTURES}")
    endif()
    set(CPACK_SYSTEM_NAME "${_charon_os}-${_charon_arch}")
endif()
set(CPACK_PACKAGE_FILE_NAME "charon-${PROJECT_VERSION}-${CPACK_SYSTEM_NAME}")
if(WIN32)
    set(CPACK_GENERATOR "ZIP")
else()
    set(CPACK_GENERATOR "TGZ")
endif()
include(CPack)
