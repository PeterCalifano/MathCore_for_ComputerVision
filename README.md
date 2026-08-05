# MathCore for Computer Vision

MathCore is a header-only C++20 mathematics library built on Eigen. It provides reusable linear-algebra, rotation, interpolation, random-sampling, and logging utilities for computer-vision and space-navigation applications.

The public CMake/package identity remains `mathcore_for_cv`. The canonical repository is [PeterCalifano/MathCore_for_ComputerVision](https://github.com/PeterCalifano/MathCore_for_ComputerVision); an older local checkout directory may still be named `MathCore_for_SpaceNav` without affecting consumers.

## Requirements

| Dependency | Version | Use |
|---|---:|---|
| CMake | 3.15+ | Configure, install, and package |
| C++ compiler | C++20 | GCC or Clang |
| Eigen | 3.4+ | Required public dependency |
| Catch2 | 3.x | Tests; fetched when unavailable locally |
| CUDA Toolkit | 12.0+ | Optional CUDA-aware builds |
| oneTBB | any | Optional parallel-runtime interface |
| OpenGL | any | Optional graphics interface |
| Python | 3.12 | Optional pytest and gtwrap bindings |
| MATLAB | R2023b+ | Optional gtwrap bindings |
| Doxygen and Graphviz | any | API documentation |

MathCore intentionally does not require ROS, OptiX, TensorRT, ZeroMQ, spdlog, or a compiled logging library.

## Quick start

```bash
# RelWithDebInfo, tests enabled
./build_lib.sh -N

# Debug or portable release builds
./build_lib.sh -N -t debug
./build_lib.sh -N -t release -D CPU_ENABLE_NATIVE_TUNING=OFF

# Explicit clean rebuild of an owned in-repository build tree
./build_lib.sh -N -B build --clean
```

The equivalent manual workflow is:

```bash
cmake -S . -B build -GNinja -DCMAKE_BUILD_TYPE=RelWithDebInfo
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

`build_lib.sh --clean` accepts only conventional in-repository `build`, `build*`, or `out/*` paths. An existing directory must contain a CMake cache owned by this checkout.

## Library layout

Public headers live below `src/mathcore/`:

- `linalg/`: vector traits, skew matrices, and validation;
- `rotations/`: quaternion traits and operations;
- `interpolation/`: interpolation primitives and Chebyshev interpolation;
- `random/`: vector sampling;
- `logging/`: the dependency-free `mathcore::logging::CLogger`;
- `mathcore.h`: umbrella include for all public modules.

The root target auto-detects whether compiled sources exist. The current MathCore layout has none, so `mathcore_for_cv` is an `INTERFACE` library and remains header-only even when logging is used.

## Important CMake options

| Option | Default | Meaning |
|---|---:|---|
| `ENABLE_TESTS` | `ON` | Register CTest tests |
| `ENABLE_PYTHON_TESTS` | `ON` | Register `test*.py` through pytest |
| `ENABLE_CUDA` | `OFF` | Enable the CUDA language and interface flags |
| `ENABLE_TBB` | `OFF` | Enable oneTBB |
| `ENABLE_OPENGL` | `OFF` | Enable OpenGL |
| `mathcore_for_cv_BUILD_PROGRAMS` | `ON` | Build `src/bin` when top-level |
| `mathcore_for_cv_BUILD_EXAMPLES` | `ON` | Build examples when top-level |
| `CPU_ENABLE_NATIVE_TUNING` | `ON` | Add native CPU tuning in optimized native builds |
| `CPU_ENABLE_SIMD` | `OFF` | Select an explicit SIMD ISA |
| `CPU_ENABLE_FMA` | `OFF` | Enable CPU FMA instructions |
| `CUDA_ENABLE_FMAD` | `ON` | Control NVCC FMA contraction |
| `CUDA_ENABLE_EXTRA_DEVICE_VECTORIZATION` | `OFF` | Enable NVCC device vectorization |
| `CUDA_USE_FAST_MATH` | `OFF` | Enable regular CUDA fast math |
| `NO_OPTIMIZATION` | `OFF` | Use profiler-friendly flags while retaining assertions |
| `WRITE_SOURCE_VERSION_FILE` | `OFF` | Opt in to writing `VERSION` in the source tree |

Cross-compilation disables host-native tuning automatically. A cross build must use target-specific flags through `CPU_EXTRA_OPT_FLAGS`; `CPU_SIMD_LEVEL=native` is rejected while cross-compiling.

## Header-only logger

```cpp
#include <mathcore/logging/logger.h>

mathcore::logging::CLogger logger{"triangulation",
    mathcore::logging::ELogLevel::Info};
logger.info("tracks: ", 42);
```

Messages at critical/error/warning severity use the diagnostic stream; info/debug/trace use the output stream. Color is disabled by default. Call `setLevelFromEnvironment()` to explicitly read `MATHCORE_FOR_CV_LOG_LEVEL`; accepted values are `quiet`, `critical`, `error`, `warning`, `info`, `debug`, `trace`, their documented aliases, or integers `0` through `6`.

## Install and consume

```bash
cmake -S . -B build-install -GNinja \
  -DCMAKE_INSTALL_PREFIX="$PWD/install" \
  -DCPU_ENABLE_NATIVE_TUNING=OFF
cmake --build build-install
cmake --install build-install
```

Downstream CMake:

```cmake
find_package(mathcore_for_cv REQUIRED)
target_link_libraries(my_target PRIVATE mathcore_for_cv::mathcore_for_cv)
```

The installed package exports Eigen and any enabled CUDA, TBB, or OpenGL dependency. Build-tree and installed `VERSION` metadata use the version resolved from Git tags, with the source-tree fallback used only when Git metadata is unavailable.

Create binary and source archives with:

```bash
cpack --config build-install/CPackConfig.cmake
cpack --config build-install/CPackSourceConfig.cmake
```

Both archives receive the authoritative build-generated `VERSION`; nested build trees and generated Python caches are excluded from source packages.

## Python and MATLAB wrappers

Wrappers use gtwrap from an explicit checkout, `lib/wrap`, `wrap`, or an installed `gtwrap` package:

```bash
# Python wrapper
./build_lib.sh -N -p

# MATLAB wrapper
MATLAB_ROOT_DIR=/usr/local/MATLAB/R2024b ./build_lib.sh -N -m

# Explicit wrapper checkout
./build_lib.sh -N -p --gtwrap-root /path/to/wrap
```

Ordinary configuration never updates, initializes, or creates the wrapper checkout. `--wrap-update` and `--wrap-submodule-init` are explicit maintenance operations. Direct CMake callers must additionally grant `GTWRAP_MAINTENANCE_UPDATE=ON` before requesting synchronization.

Python `setup.py`, `pyproject.toml`, the package copy, and `_wrapper_build.py` are generated only below the CMake build tree. Build a wheel from that configured directory:

```bash
python -m pip wheel build/python --no-build-isolation --no-deps
```

The wheel stages only the exact wrapper and project-owned shared runtimes declared by CMake, using loader-relative runtime paths. Checkout-only `_wrapper_build.py` metadata is not installed.

MATLAB library compatibility can be inspected without mutation:

```bash
./scripts/use_system_matlab_libraries.sh --matlab-version R2024b --all
sudo ./scripts/use_system_matlab_libraries.sh --matlab-version R2024b --all --apply
sudo ./scripts/use_system_matlab_libraries.sh --matlab-version R2024b --all --restore
```

The utility never invokes `sudo` itself. Mutation requires an explicit mode and root privileges, prints every command/output/status, and retains one-time recovery symlinks.

## Documentation

```bash
cmake --preset docs
cmake --build --preset docs
```

HTML and XML are written below `build_docs/doc/`. `DOC_WARN_AS_ERROR=ON` is available for strict documentation gates; the default allows legacy documentation warnings while still producing the site.

## Development container

The checked-in devcontainer uses Ubuntu 24.04, Python 3.12, and CUDA 12.9. Regenerate configuration without building or running a container:

```bash
./configure_devcontainer.sh --cuda --gpu-runtime docker
./configure_devcontainer.sh --no-cuda --base ubuntu-24.04
```

Docker and Podman GPU argument styles are supported. Standalone Dockerfile builds can request the CUDA apt setup with `--build-arg INSTALL_CUDA=on`; this is separate from the normal devcontainer feature path.

## Tests and CI

CTest is the common entry point. Catch2 cases receive the `catch2` label, and pytest files receive `python;pytest` labels.

```bash
ctest --test-dir build --output-on-failure
ctest --test-dir build --output-on-failure -L catch2
./tests/scripts/test_use_system_matlab_libraries.sh
```

GitHub workflows cover hosted Linux CPU builds, explicitly labeled self-hosted CUDA builds, and Doxygen/Pages generation. They use portable CPU flags; release artifacts intended for a particular host may re-enable native tuning explicitly.
