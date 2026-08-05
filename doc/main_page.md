# MathCore for Computer Vision {#mainpage}

MathCore is a header-only C++20 Eigen library for linear algebra, rotations,
interpolation, random sampling, and dependency-free component logging. The
public CMake/package identity is `mathcore_for_cv`.

## Build and test

```bash
./build_lib.sh -N
ctest --test-dir build --output-on-failure
```

## Optional capabilities

```bash
# CUDA-aware configuration
./build_lib.sh -N -D ENABLE_CUDA=ON

# oneTBB and explicit SIMD/FMA
./build_lib.sh -N -D ENABLE_TBB=ON \
  -D CPU_ENABLE_SIMD=ON -D CPU_SIMD_LEVEL=avx2 -D CPU_ENABLE_FMA=ON

# Python and MATLAB wrappers
./build_lib.sh -N -p
MATLAB_ROOT_DIR=/usr/local/MATLAB/R2024b ./build_lib.sh -N -m
```

## Downstream CMake

```cmake
find_package(mathcore_for_cv REQUIRED)
target_link_libraries(my_target PRIVATE mathcore_for_cv::mathcore_for_cv)
```

See the repository `README.md` for versioning, packaging, wrapper-maintenance,
logger, and devcontainer contracts.
