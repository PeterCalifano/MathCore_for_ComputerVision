# CLAUDE.md

## Project overview

MathCore for Computer Vision is a header-only C++20 Eigen library. Its stable CMake/package identity is `mathcore_for_cv`, regardless of the checkout directory name. Public modules cover linear algebra, rotations, interpolation, random sampling, and `mathcore::logging::CLogger`.

## Common commands

```bash
./build_lib.sh -N
./build_lib.sh -N -t debug
./build_lib.sh -N -t release -D CPU_ENABLE_NATIVE_TUNING=OFF
ctest --test-dir build --output-on-failure

cmake --preset docs
cmake --build --preset docs
```

Optional CUDA, Python wrapper, and MATLAB wrapper builds:

```bash
./build_lib.sh -N -D ENABLE_CUDA=ON
./build_lib.sh -N -p
MATLAB_ROOT_DIR=/usr/local/MATLAB/R2024b ./build_lib.sh -N -m
```

## Architecture and ownership

- Public headers are below `src/mathcore/`; `src/mathcore/mathcore.h` is the umbrella header.
- Each module owns a `CMakeLists.txt` that installs its public headers. Do not add a `.cpp` merely for convenience: any compiled source changes the root target from `INTERFACE` to a binary library.
- Tests named `tests/test*.cpp` are discovered as Catch2 executables. Python `test*.py` files are registered through pytest when enabled.
- Generic build behavior belongs in `cmake/`; product-specific dependency decisions stay in the root and `src/CMakeLists.txt`.
- `lib/wrap` is an externally owned gitlink. Ordinary configuration must not move it.

## Wrapper and packaging contracts

- Wrapper checkout updates and submodule initialization are explicit maintenance operations only.
- Python `setup.py`, `pyproject.toml`, package copies, and `_wrapper_build.py` are generated below `<build>/python`, never in the source tree.
- `_wrapper_build.py` is checkout-only metadata and must not enter installs or wheels.
- Python installs remain relative to `CMAKE_INSTALL_PREFIX`; pip owns active-environment installation.
- The MATLAB host-library utility is dry-run by default and never invokes sudo internally.

## Scope boundaries

MathCore does not carry ROS, OptiX/PTX, TensorRT, ZeroMQ, or spdlog facilities. Keep CUDA optional and OFF by default. Container changes require static configuration validation; image/runtime testing is not part of the normal repository gate.
