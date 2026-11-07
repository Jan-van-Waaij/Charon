#pragma once

#include <cstdint>
#include <mutex>
#include <string>

namespace charon {

/// The shortest decimal representation that reads back as the same double, in the
/// style of Julia's print: "0.25", "1.0", "1.0e-5", "-Inf", "NaN".
std::string format_number(double x);

/// Output file; written gzip-compressed when the path ends in ".gz".
class OutputFile {
public:
    explicit OutputFile(const std::string& path);
    ~OutputFile();
    OutputFile(const OutputFile&) = delete;
    OutputFile& operator=(const OutputFile&) = delete;

    void write(const std::string& text);
    /// Flush and close; throws on write errors.
    void close();

private:
    void flush_buffer();

    std::string path_;
    void* file_ = nullptr;  // gzFile
    std::string buffer_;
};

/// Thread-safe messages on standard error.
class Logger {
public:
    explicit Logger(bool quiet = false) : quiet_(quiet) {}
    void info(const std::string& message);
    /// Messages that are printed even in quiet mode.
    void warn(const std::string& message);
    bool quiet() const { return quiet_; }

private:
    bool quiet_;
    std::mutex mutex_;
};

}  // namespace charon
