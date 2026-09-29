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
./build_lib.sh -N -D mathcore_for_cv_ENABLE_CUDA=ON

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

## MATLAB utilities

Load MATLAB source helpers directly through the MATLAB path. They do not require
generated C++ bindings.

### File and data hashing

Use Jan Simon's bundled DataHash package for file and MATLAB data hashes.
The package was imported from the installed add-on, version 1.7.1, on
29 September 2026: `DataHash.m`, `uTest_DataHash.m` and the
[BSD license](../matlab/misc/DataHash/license.txt). The source headers identify revision 043
(18 April 2019) and test revision 004 (19 May 2019).

Keep `DataHash.m` and `license.txt` byte-for-byte unchanged. Retain the upstream
documentation and formatting in the test; change its file-fixture write to
`fwrite(fid, uint8(Test{1}), 'uint8')` so character values above 127 remain raw
bytes under MATLAB's UTF-8 default. The original all-byte fixture wrote 384 bytes
while its expected digest described 256 bytes. Record this compatibility fix
in the test header. The local `.gitattributes` exempts only the two vendored
MATLAB files from end-of-line whitespace checks.
Exclude the add-on's `resources/` directory, which contains installer metadata,
an archive and screenshots. The copied package source was
`/home/peterc/MATLAB-Add-Ons/Collections/DataHash`; that path is provenance only
and is never used at runtime.

Use `matlab/misc/ComputeFileSha256.m` for raw-file SHA-256.
It preserves the existing caller contract and invokes
`DataHash(path, 'SHA-256', 'file', 'hex')`. File mode reads bounded chunks and
returns a lowercase hexadecimal digest. The default `DataHash(data)` call
instead hashes MATLAB array contents, type and dimensions using MD5; select
the input mode explicitly when verifying files. Both calls require the MATLAB
JVM and run on the host, outside MATLAB Coder paths.

Run this example from the MathCore repository:

```matlab
addpath(genpath('matlab/misc'));
charFilePath = tempname;
dFileId = fopen(charFilePath, 'wb');
fwrite(dFileId, uint8('abc'), 'uint8');
fclose(dFileId);
fprintf('%s\n', ComputeFileSha256(charFilePath));
delete(charFilePath);
uTest_DataHash(false);
```

The first output line is:

```text
ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad
```

The upstream test prints its checks and raises an error on failure. Pass `false`
to skip its speed benchmark.
