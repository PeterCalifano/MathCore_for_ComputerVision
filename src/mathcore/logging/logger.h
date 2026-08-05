#pragma once

/// @file
/// @brief Dependency-free, header-only logging for MathCore components.

#include <algorithm>
#include <atomic>
#include <cctype>
#include <concepts>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <mutex>
#include <optional>
#include <ostream>
#include <sstream>
#include <string>
#include <string_view>
#include <utility>

namespace mathcore::logging
{
    /** @brief Ordered logging threshold used by CLogger. */
    enum class ELogLevel : std::uint8_t
    {
        Quiet = 0,
        Critical = 1,
        Error = 2,
        Warning = 3,
        Info = 4,
        Debug = 5,
        Trace = 6
    };

    /** @brief Select whether CLogger emits ANSI terminal color sequences. */
    enum class ELogColorMode : std::uint8_t
    {
        Disabled = 0,
        Enabled = 1
    };

    /** @brief A value that can be appended to a standard output stream. */
    template <typename TValue>
    concept StreamInsertable = requires(std::ostream &stream, TValue &&value) {
        stream << std::forward<TValue>(value);
    };

    /**
     * @brief Small dependency-free, component-scoped logger.
     *
     * Messages are assembled before acquiring a process-wide output lock, so
     * complete lines from independent logger instances cannot interleave.
     * Diagnostic severities use the diagnostic stream; informational
     * severities use the ordinary output stream.
     */
    class CLogger final
    {
      public:
        /**
         * @brief Construct a logger for one component.
         *
         * @param componentName Component printed in each line. Empty selects
         * `mathcore_for_cv`.
         * @param level Initial verbosity threshold.
         * @param colorMode Explicit ANSI color policy.
         * @param outputStream Stream for info, debug, and trace messages.
         * @param diagnosticStream Stream for critical, error, and warning messages.
         */
        explicit CLogger(std::string componentName = "mathcore_for_cv",
                         ELogLevel level = ELogLevel::Error,
                         ELogColorMode colorMode = ELogColorMode::Disabled,
                         std::ostream &outputStream = std::cout,
                         std::ostream &diagnosticStream = std::clog)
            : componentName_(componentName.empty() ? "mathcore_for_cv"
                                                   : std::move(componentName)),
              level_(level), colorMode_(colorMode), outputStream_(outputStream),
              diagnosticStream_(diagnosticStream)
        {
        }

        CLogger(const CLogger &) = delete;
        CLogger &operator=(const CLogger &) = delete;
        CLogger(CLogger &&) = delete;
        CLogger &operator=(CLogger &&) = delete;
        ~CLogger() = default;

        /** @brief Change the active verbosity threshold. */
        void setLevel(const ELogLevel level) noexcept
        {
            level_.store(level, std::memory_order_relaxed);
        }

        /** @brief Return the active verbosity threshold. */
        [[nodiscard]] ELogLevel getLevel() const noexcept
        {
            return level_.load(std::memory_order_relaxed);
        }

        /** @brief Return true when a severity is enabled by the active threshold. */
        [[nodiscard]] bool shouldLog(const ELogLevel severity) const noexcept
        {
            const ELogLevel configuredLevel = getLevel();
            const auto configuredValue = static_cast<std::uint8_t>(configuredLevel);
            const auto severityValue = static_cast<std::uint8_t>(severity);

            if (configuredLevel == ELogLevel::Quiet || severity == ELogLevel::Quiet)
            {
                return false;
            }
            if (configuredValue > static_cast<std::uint8_t>(ELogLevel::Trace) ||
                severityValue > static_cast<std::uint8_t>(ELogLevel::Trace))
            {
                return false;
            }
            return severityValue <= configuredValue;
        }

        /**
         * @brief Parse a case-insensitive level name or numeric value from 0 to 6.
         *
         * Leading and trailing ASCII whitespace is ignored. Accepted aliases are
         * quiet/off, critical/fatal, warning/warn, error, info, debug, and trace.
         */
        [[nodiscard]] static std::optional<ELogLevel> tryParseLevel(std::string_view text)
        {
            text = trimAsciiWhitespace_(text);
            if (text.size() == 1 && text.front() >= '0' && text.front() <= '6')
            {
                return static_cast<ELogLevel>(text.front() - '0');
            }

            std::string normalized(text);
            std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                           [](const unsigned char value)
                           {
                               return static_cast<char>(std::tolower(value));
                           });

            if (normalized == "quiet" || normalized == "off")
            {
                return ELogLevel::Quiet;
            }
            if (normalized == "critical" || normalized == "fatal")
            {
                return ELogLevel::Critical;
            }
            if (normalized == "error")
            {
                return ELogLevel::Error;
            }
            if (normalized == "warning" || normalized == "warn")
            {
                return ELogLevel::Warning;
            }
            if (normalized == "info")
            {
                return ELogLevel::Info;
            }
            if (normalized == "debug")
            {
                return ELogLevel::Debug;
            }
            if (normalized == "trace")
            {
                return ELogLevel::Trace;
            }
            return std::nullopt;
        }

        /**
         * @brief Apply a valid level from an environment variable.
         * @return True only when a valid value was found and applied.
         */
        bool setLevelFromEnvironment(
            const std::string_view variableName = "MATHCORE_FOR_CV_LOG_LEVEL")
        {
            if (variableName.empty())
            {
                return false;
            }

            const std::string variableNameCopy(variableName);
            const char *environmentValue = std::getenv(variableNameCopy.c_str());
            if (environmentValue == nullptr)
            {
                return false;
            }

            const auto parsedLevel = tryParseLevel(environmentValue);
            if (!parsedLevel.has_value())
            {
                return false;
            }
            setLevel(*parsedLevel);
            return true;
        }

        /** @brief Emit a critical message when enabled. */
        template <StreamInsertable... TArgs>
        void critical(TArgs &&...args)
        {
            write_(ELogLevel::Critical, std::forward<TArgs>(args)...);
        }

        /** @brief Emit an error message when enabled. */
        template <StreamInsertable... TArgs>
        void error(TArgs &&...args)
        {
            write_(ELogLevel::Error, std::forward<TArgs>(args)...);
        }

        /** @brief Emit a warning message when enabled. */
        template <StreamInsertable... TArgs>
        void warning(TArgs &&...args)
        {
            write_(ELogLevel::Warning, std::forward<TArgs>(args)...);
        }

        /** @brief Emit an informational message when enabled. */
        template <StreamInsertable... TArgs>
        void info(TArgs &&...args)
        {
            write_(ELogLevel::Info, std::forward<TArgs>(args)...);
        }

        /** @brief Emit a debug message when enabled. */
        template <StreamInsertable... TArgs>
        void debug(TArgs &&...args)
        {
            write_(ELogLevel::Debug, std::forward<TArgs>(args)...);
        }

        /** @brief Emit a trace message when enabled. */
        template <StreamInsertable... TArgs>
        void trace(TArgs &&...args)
        {
            write_(ELogLevel::Trace, std::forward<TArgs>(args)...);
        }

      private:
        template <StreamInsertable... TArgs>
        void write_(const ELogLevel severity, TArgs &&...args)
        {
            if (!shouldLog(severity))
            {
                return;
            }

            std::ostringstream messageStream;
            (messageStream << ... << std::forward<TArgs>(args));
            writeMessage_(severity, messageStream.str());
        }

        void writeMessage_(const ELogLevel severity, const std::string_view message)
        {
            std::ostringstream line;
            if (colorMode_ == ELogColorMode::Enabled)
            {
                line << colorCode_(severity);
            }
            line << '[' << componentName_ << "][" << levelLabel_(severity) << "] " << message;
            if (colorMode_ == ELogColorMode::Enabled)
            {
                line << "\033[0m";
            }
            line << '\n';

            const std::scoped_lock outputLock(outputMutex_());
            selectStream_(severity) << line.str();
        }

        [[nodiscard]] std::ostream &selectStream_(const ELogLevel severity) const noexcept
        {
            if (severity == ELogLevel::Critical || severity == ELogLevel::Error ||
                severity == ELogLevel::Warning)
            {
                return diagnosticStream_;
            }
            return outputStream_;
        }

        [[nodiscard]] static std::string_view trimAsciiWhitespace_(std::string_view text)
        {
            while (!text.empty() &&
                   std::isspace(static_cast<unsigned char>(text.front())) != 0)
            {
                text.remove_prefix(1);
            }
            while (!text.empty() &&
                   std::isspace(static_cast<unsigned char>(text.back())) != 0)
            {
                text.remove_suffix(1);
            }
            return text;
        }

        [[nodiscard]] static std::mutex &outputMutex_()
        {
            static std::mutex outputMutex;
            return outputMutex;
        }

        [[nodiscard]] static std::string_view levelLabel_(const ELogLevel severity) noexcept
        {
            switch (severity)
            {
            case ELogLevel::Critical:
                return "CRITICAL";
            case ELogLevel::Error:
                return "ERROR";
            case ELogLevel::Warning:
                return "WARNING";
            case ELogLevel::Info:
                return "INFO";
            case ELogLevel::Debug:
                return "DEBUG";
            case ELogLevel::Trace:
                return "TRACE";
            case ELogLevel::Quiet:
            default:
                return "QUIET";
            }
        }

        [[nodiscard]] static std::string_view colorCode_(const ELogLevel severity) noexcept
        {
            switch (severity)
            {
            case ELogLevel::Critical:
                return "\033[1;31m";
            case ELogLevel::Error:
                return "\033[31m";
            case ELogLevel::Warning:
                return "\033[33m";
            case ELogLevel::Info:
                return "\033[34m";
            case ELogLevel::Debug:
                return "\033[36m";
            case ELogLevel::Trace:
                return "\033[2m";
            case ELogLevel::Quiet:
            default:
                return {};
            }
        }

      private:
        std::string componentName_;
        std::atomic<ELogLevel> level_;
        ELogColorMode colorMode_;
        std::ostream &outputStream_;
        std::ostream &diagnosticStream_;
    };
} // namespace mathcore::logging
