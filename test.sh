#!/usr/bin/env bash
# Builds and runs every test suite in tests/.   Usage: ./test.sh
# Set BUILD_DIR to put the test binaries somewhere other than ./build.
cd "$(dirname "$0")" || exit 1

BUILD_DIR="${BUILD_DIR:-build}"
CXX="${CXX:-g++}"
CC="${CC:-gcc}"
FLAGS="-std=c++20 -Wall -Wextra -Wpedantic -Iinclude"
LIB_SRC=$(ls src/*.cpp)          # all library sources, so new files are picked up automatically

mkdir -p "$BUILD_DIR"
status=0

# Shewchuk's predicates are C: compile once to an object file, then link it.
$CC -w -c src/predicate.c -o "$BUILD_DIR/predicate.o" || exit 1

for t in test_geometry test_delaunay test_io test_domain test_refine; do
    echo "===== $t ====="
    if $CXX $FLAGS tests/$t.cpp $LIB_SRC "$BUILD_DIR/predicate.o" -o "$BUILD_DIR/$t"; then
        "$BUILD_DIR/$t" || status=1
    else
        echo "BUILD FAILED: $t"
        status=1
    fi
    echo
done

if [ $status -eq 0 ]; then echo "ALL TEST SUITES PASSED"; else echo "SOME TESTS FAILED"; fi
exit $status
