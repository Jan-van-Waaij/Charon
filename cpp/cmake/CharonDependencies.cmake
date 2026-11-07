# Provides the imported targets charon::eigen (header only) and charon::zlib.
#
# By default the system packages are used (Debian/Ubuntu: libeigen3-dev zlib1g-dev,
# Fedora: eigen3-devel zlib-devel / zlib-static, macOS: brew install eigen,
# MSYS2: mingw-w64-x86_64-eigen3 mingw-w64-x86_64-zlib).
# With -DCHARON_FETCH_DEPS=ON both are downloaded and zlib is built from source,
# which only needs a C compiler and network access at configure time.

include(FetchContent)
# Keep FetchContent_Populate (only the sources are needed, not the projects) without
# deprecation warnings in newer CMake, and use the archive timestamps policy of 3.24+.
if(POLICY CMP0169)
    cmake_policy(SET CMP0169 OLD)
endif()
if(POLICY CMP0135)
    cmake_policy(SET CMP0135 NEW)
endif()

set(CHARON_EIGEN_URL "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz")
set(CHARON_EIGEN_SHA256 "8586084f71f9bde545ee7fa6d00288b264a2b7ac3607b974e54d13e7162c1c72")
set(CHARON_ZLIB_URL "https://github.com/madler/zlib/releases/download/v1.3.1/zlib-1.3.1.tar.gz")
set(CHARON_ZLIB_SHA256 "9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23")

# ---------------------------------------------------------------- Eigen
if(NOT CHARON_FETCH_DEPS)
    find_package(Eigen3 3.3 QUIET NO_MODULE)
endif()
if(TARGET Eigen3::Eigen)
    add_library(charon_eigen INTERFACE)
    target_link_libraries(charon_eigen INTERFACE Eigen3::Eigen)
    message(STATUS "charon: using system Eigen ${Eigen3_VERSION}")
else()
    if(NOT CHARON_FETCH_DEPS)
        message(STATUS "charon: Eigen not found on the system, downloading Eigen 3.4.0")
    endif()
    FetchContent_Declare(charon_eigen_src URL "${CHARON_EIGEN_URL}" URL_HASH SHA256=${CHARON_EIGEN_SHA256})
    FetchContent_GetProperties(charon_eigen_src)
    if(NOT charon_eigen_src_POPULATED)
        # Only the headers are needed; do not add Eigen's own CMake project.
        FetchContent_Populate(charon_eigen_src)
    endif()
    add_library(charon_eigen INTERFACE)
    target_include_directories(charon_eigen SYSTEM INTERFACE "${charon_eigen_src_SOURCE_DIR}")
endif()
add_library(charon::eigen ALIAS charon_eigen)

# ---------------------------------------------------------------- zlib
if(NOT CHARON_FETCH_DEPS)
    if(CHARON_STATIC AND NOT MSVC AND NOT ZLIB_LIBRARY)
        # Prefer the static archive, a fully static link cannot use libz.so.
        find_library(ZLIB_LIBRARY NAMES libz.a z.a libzlibstatic.a zlibstatic)
    endif()
    find_package(ZLIB QUIET)
endif()
if(ZLIB_FOUND)
    add_library(charon_zlib INTERFACE)
    target_link_libraries(charon_zlib INTERFACE ZLIB::ZLIB)
    message(STATUS "charon: using system zlib ${ZLIB_VERSION_STRING} (${ZLIB_LIBRARIES})")
else()
    if(NOT CHARON_FETCH_DEPS)
        message(STATUS "charon: zlib not found on the system, building zlib 1.3.1 from source")
    endif()
    enable_language(C)
    FetchContent_Declare(charon_zlib_src URL "${CHARON_ZLIB_URL}" URL_HASH SHA256=${CHARON_ZLIB_SHA256})
    FetchContent_GetProperties(charon_zlib_src)
    if(NOT charon_zlib_src_POPULATED)
        FetchContent_Populate(charon_zlib_src)
    endif()
    # Compile the zlib sources directly into a static library; this avoids zlib's
    # own CMake project (which also builds a shared library and example programs).
    set(_zsrc "${charon_zlib_src_SOURCE_DIR}")
    add_library(charon_zlib STATIC
        ${_zsrc}/adler32.c ${_zsrc}/compress.c ${_zsrc}/crc32.c ${_zsrc}/deflate.c
        ${_zsrc}/gzclose.c ${_zsrc}/gzlib.c ${_zsrc}/gzread.c ${_zsrc}/gzwrite.c
        ${_zsrc}/infback.c ${_zsrc}/inffast.c ${_zsrc}/inflate.c ${_zsrc}/inftrees.c
        ${_zsrc}/trees.c ${_zsrc}/uncompr.c ${_zsrc}/zutil.c)
    target_include_directories(charon_zlib SYSTEM PUBLIC "${_zsrc}")
    if(NOT WIN32)
        target_compile_definitions(charon_zlib PRIVATE HAVE_UNISTD_H)
    endif()
    set_target_properties(charon_zlib PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()
add_library(charon::zlib ALIAS charon_zlib)
