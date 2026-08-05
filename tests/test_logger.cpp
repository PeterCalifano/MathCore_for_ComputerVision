#include <catch2/catch_test_macros.hpp>

#include <mathcore/logging/logger.h>

#include <cstdlib>
#include <sstream>
#include <string>

using mathcore::logging::CLogger;
using mathcore::logging::ELogColorMode;
using mathcore::logging::ELogLevel;

TEST_CASE("logger_routes_enabled_severities_to_the_expected_streams", "[logging]")
{
    std::ostringstream output;
    std::ostringstream diagnostic;
    CLogger logger{"solver", ELogLevel::Info, ELogColorMode::Disabled, output, diagnostic};

    logger.error("failed ", 7);
    logger.info("ready");
    logger.debug("hidden");

    REQUIRE(diagnostic.str() == "[solver][ERROR] failed 7\n");
    REQUIRE(output.str() == "[solver][INFO] ready\n");
}

TEST_CASE("logger_parses_named_and_numeric_levels", "[logging]")
{
    REQUIRE(CLogger::tryParseLevel(" Warn ") == ELogLevel::Warning);
    REQUIRE(CLogger::tryParseLevel("6") == ELogLevel::Trace);
    REQUIRE_FALSE(CLogger::tryParseLevel("verbose").has_value());
    REQUIRE_FALSE(CLogger{}.shouldLog(ELogLevel::Info));
}

TEST_CASE("logger_uses_mathcore_defaults", "[logging]")
{
    std::ostringstream output;
    std::ostringstream diagnostic;
    CLogger logger{"", ELogLevel::Info, ELogColorMode::Disabled, output, diagnostic};

    logger.info("ready");
    REQUIRE(output.str() == "[mathcore_for_cv][INFO] ready\n");

    REQUIRE(::setenv("MATHCORE_FOR_CV_LOG_LEVEL", "trace", 1) == 0);
    REQUIRE(logger.setLevelFromEnvironment());
    REQUIRE(logger.getLevel() == ELogLevel::Trace);
    REQUIRE(::unsetenv("MATHCORE_FOR_CV_LOG_LEVEL") == 0);
}
