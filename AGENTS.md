# Repository operating guide

MathCore is a header-only C++20 library with stable public identity `mathcore_for_cv`. Preserve that identity in targets, exported packages, wrappers, and include paths; the human-facing repository name is MathCore for Computer Vision.

## Implementation rules

- Keep public code below `src/mathcore/` and in namespace `mathcore`.
- Use Doxygen file and public API documentation for C++ headers.
- Prefer concepts to SFINAE and Catch2 for C++ behavior tests.
- Add new public headers to the owning module install rules and, when appropriate, `mathcore/mathcore.h`.
- Do not introduce compiled sources for header-only utilities.

## Build and verification

```bash
./build_lib.sh -N
ctest --test-dir build --output-on-failure
cmake --preset docs
cmake --build --preset docs
```

Use fresh out-of-tree builds for CMake options, install/export, packaging, cross-toolchain, and wrapper acceptance. Do not import donor `VerifyTemplateProject*` self-conformance tests into this product suite.

## Safety boundaries

- Treat `lib/wrap` as external: no checkout update, submodule initialization, or gitlink change during ordinary configure/build work.
- Generated wrapper metadata belongs below the build tree, not in `python/` sources.
- `build_lib.sh --clean` may remove only an owned, conventional in-repository build directory.
- Do not add ROS, OptiX/PTX, TensorRT, ZeroMQ, or spdlog facilities without an explicit product decision.
- Container validation is configuration-only unless runtime testing is explicitly requested.

## Review handoff

Before presenting a staged batch, inspect `git diff --cached`, confirm new/modified public APIs have file-level and callable documentation, and keep unrelated user work outside the index. Do not commit or push unless asked.
