# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MathCore_for_SpaceNav is a **header-only C++20 math library** built on Eigen, providing reusable utilities for linear algebra, rotations, interpolation, and random sampling. The CMake project name is `mathcore_for_cv`.

## Build Commands

```bash
# Full configure + build + test (default: RelWithDebInfo)
./build_lib.sh

# Debug build
./build_lib.sh -t debug

# Release build (forces tests)
./build_lib.sh -t release

# Rebuild only (skip CMake configure)
./build_lib.sh -r

# Build with Ninja, 8 jobs
./build_lib.sh -t debug -j 8 -N

# Build and install
./build_lib.sh -i

# Clean rebuild
./build_lib.sh --clean
```

Manual CMake workflow:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_EXPORT_COMPILE_COMMANDS=ON
cmake --build build --parallel $(nproc)
```

## Running Tests

Tests use **Catch2** (fetched automatically if not found).

```bash
# Run all tests
ctest --test-dir build --output-on-failure

# Run a single test binary
./build/tests/test_linalg

# Run specific Catch2 test cases by name
./build/tests/test_linalg "[linalg]"
```

Test files follow the naming pattern `tests/test_<module>.cpp`. Each test file links against `Catch2::Catch2WithMain`. The `add_tests()` CMake function in `cmake/cmake_utils.cmake` auto-discovers test files in `tests/`.

## Architecture

### Library structure (all under `src/mathcore/`)

- **`mathcore.h`** — Umbrella header that includes all public modules
- **`linalg/`** — Skew-symmetric matrices, SPD validation, vector traits
- **`rotations/`** — Quaternion operations and traits (Eigen-based)
- **`interpolation/`** — Polynomial interpolation base class, Chebyshev interpolation
- **`random/`** — Random scatter/sampling utilities

All code lives in the `mathcore` namespace. Headers use `#pragma once` guards.

### Key dependencies

- **Eigen 3.4+** (required) — Core linear algebra types used throughout
- **Catch2** — Testing framework (auto-fetched via `cmake/HandleCatch2.cmake`)
- **C++20** — Uses concepts, `<concepts>`, and modern standard library features

### CMake design patterns

- The library auto-detects whether it's header-only or compiled based on presence of `.cpp` files in `src/`. Currently header-only (INTERFACE target).
- Interface targets with `${LIB_NAMESPACE}_` prefix are used for compile options (CUDA, OpenGL, sanitizers, profiling, etc.) to avoid global target name clashes.
- Submodules in `lib/` are auto-discovered by `cmake/HandleSubmodules.cmake`.
- The `lib/wrap/` submodule (gtwrap + pybind11) enables optional Python/MATLAB wrapper generation via `--python-wrap` / `--matlab-wrap` build flags.
- Cross-compilation toolchains are in `cmake/toolchains/defaults/`.

### Extensibility pattern

New modules go in `src/mathcore/<module>/` with their own `CMakeLists.txt` that installs headers and appends sources to the parent scope. Add the subdirectory in `src/mathcore/CMakeLists.txt` and include the new header in the umbrella `src/mathcore/mathcore.h`.
