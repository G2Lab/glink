cmake_minimum_required(VERSION 3.18)
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "Build outside the source tree: cmake -B build")
endif()
project(glink3 LANGUAGES C CXX)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Stores are HDF5 files, read and written through HighFive. Both are fetched at pinned versions and linked statically.
include(FetchContent)
foreach(off BUILD_SHARED_LIBS BUILD_TESTING HDF5_BUILD_TOOLS HDF5_BUILD_HL_LIB HDF5_BUILD_EXAMPLES HDF5_BUILD_CPP_LIB
        HDF5_ENABLE_Z_LIB_SUPPORT HDF5_ENABLE_SZIP_SUPPORT HIGHFIVE_FIND_HDF5 HIGHFIVE_UNIT_TESTS HIGHFIVE_EXAMPLES HIGHFIVE_BUILD_DOCS)
    set(${off} OFF CACHE BOOL "" FORCE)
endforeach()
FetchContent_Declare(hdf5
    URL https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5_1.14.6.tar.gz
    URL_HASH SHA256=09ee1c671a87401a5201c06106650f62badeea5a3b3941e9b1e2e1e08317357f
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_Declare(highfive
    URL https://github.com/highfive-devs/highfive/archive/refs/tags/v3.1.1.tar.gz
    URL_HASH SHA256=622034f34badda41255d7793e1c5a3046954dcf0875b0bca076e7c77088a8890
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE)
FetchContent_MakeAvailable(hdf5 highfive)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")  # no -DNDEBUG
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
add_compile_options(-Wall -Wno-sign-compare)

# Leaks and undefined behaviour are errors; floating point results must not change.
option(GLINK_SANITIZE "Build under AddressSanitizer and UBSan" OFF)
if(GLINK_SANITIZE)
    add_compile_options(-O1 -g -fno-omit-frame-pointer -fsanitize=address,undefined -fno-sanitize-recover=all)
    add_link_options(-fsanitize=address,undefined)
endif()

find_package(ZLIB REQUIRED)
find_package(OpenSSL REQUIRED)  # SHA-256, with hardware acceleration where the CPU has it
find_package(OpenMP)

add_library(glinkcore STATIC
    src/util.cpp src/store.cpp src/gfa.cpp src/coordinates.cpp
    src/conditioning.cpp src/matrix.cpp src/ld.cpp src/blocks.cpp src/matching.cpp src/index.cpp
    src/gabriel.cpp src/dosage.cpp)
target_include_directories(glinkcore PUBLIC include)
target_include_directories(glinkcore SYSTEM PUBLIC third_party)  # nlohmann/json 3.12.0
target_compile_definitions(glinkcore PUBLIC JSON_USE_IMPLICIT_CONVERSIONS=0)
target_include_directories(glinkcore SYSTEM PUBLIC ${highfive_SOURCE_DIR}/include)
target_link_libraries(glinkcore PUBLIC hdf5-static)
target_link_libraries(glinkcore PUBLIC ZLIB::ZLIB)
target_link_libraries(glinkcore PUBLIC OpenSSL::Crypto)
if(OpenMP_CXX_FOUND)
    target_link_libraries(glinkcore PUBLIC OpenMP::OpenMP_CXX)
endif()
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-mpopcnt GLINK_HAS_POPCNT)
if(GLINK_HAS_POPCNT)
    target_compile_options(glinkcore PUBLIC -mpopcnt)
endif()

add_executable(glink src/main_glink.cpp)
target_link_libraries(glink glinkcore)
add_executable(glink-compare-blocks src/main_compare_blocks.cpp)
target_link_libraries(glink-compare-blocks glinkcore)

# vg finds snarls; a run that reuses an index needs none, but building one does.
find_program(VG_EXECUTABLE NAMES vg)
if(VG_EXECUTABLE)
    message(STATUS "Found vg at: ${VG_EXECUTABLE}")
else()
    message(FATAL_ERROR "vg is required to run GLINK")
endif()

# Compares GlinkV3 with what GlinkV2 output, recorded in tests/golden by tests/make_golden.py.
enable_testing()
set(GTEST tests/third_party/googletest)
add_executable(glink-tests tests/cpp/test_graph.cpp tests/cpp/test_coordinates.cpp tests/cpp/test_conditioning.cpp
    tests/cpp/test_matrix.cpp tests/cpp/test_ld.cpp tests/cpp/test_ld_format.cpp tests/cpp/test_blocks.cpp
    tests/cpp/test_cli.cpp tests/cpp/test_index.cpp tests/cpp/test_store.cpp
    ${GTEST}/src/gtest-all.cc ${GTEST}/src/gtest_main.cc)
set_source_files_properties(${GTEST}/src/gtest-all.cc ${GTEST}/src/gtest_main.cc PROPERTIES COMPILE_OPTIONS -w)
target_include_directories(glink-tests PRIVATE tests/cpp ${GTEST})
target_include_directories(glink-tests SYSTEM PRIVATE ${GTEST}/include)
target_compile_definitions(glink-tests PRIVATE GLINK_TEST_DIR="${CMAKE_CURRENT_SOURCE_DIR}/tests"
    GLINK_BIN_DIR="$<TARGET_FILE_DIR:glink>" GLINK_VG="${VG_EXECUTABLE}")
add_dependencies(glink-tests glink glink-compare-blocks)
find_package(Threads REQUIRED)
target_link_libraries(glink-tests glinkcore Threads::Threads)
add_test(NAME unit COMMAND glink-tests)
