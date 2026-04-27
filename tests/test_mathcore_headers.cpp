#include <catch2/catch_test_macros.hpp>

#include <mathcore/interpolation/chebyshev.h>
#include <mathcore/interpolation/interpolator.h>
#include <mathcore/linalg/skew_symmetric.h>
#include <mathcore/linalg/validation.h>
#include <mathcore/mathcore.h>
#include <mathcore/random/vector_sampling.h>
#include <mathcore/rotations/quaternion.h>
#include <mathcore/rotations/quaternion_traits.h>

TEST_CASE("public_headers_compile", "[headers]")
{
    SUCCEED();
}
