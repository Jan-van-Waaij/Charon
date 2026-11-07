#include "charon/io.hpp"

#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include <zlib.h>

namespace charon {

std::string format_number(double x) {
    if (std::isnan(x)) return "NaN";
    if (std::isinf(x)) return x > 0 ? "Inf" : "-Inf";
    if (x == 0.0) return std::signbit(x) ? "-0.0" : "0.0";
    char buf[40];
    int precision = 15;  // %g drops trailing zeros, so 15 digits already give the shortest form when it has <= 15 digits
    for (; precision <= 17; ++precision) {
        std::snprintf(buf, sizeof buf, "%.*g", precision, x);
        if (std::strtod(buf, nullptr) == x) break;
    }
    std::string s(buf);
    const std::size_t e = s.find('e');
    std::string mantissa = e == std::string::npos ? s : s.substr(0, e);
    if (mantissa.find('.') == std::string::npos) mantissa += ".0";
    if (e == std::string::npos) return mantissa;
    // "1e-05" -> "1.0e-5"
    const int exponent = std::atoi(s.c_str() + e + 1);
    return mantissa + "e" + std::to_string(exponent);
}

OutputFile::OutputFile(const std::string& path) : path_(path) {
    const bool gz = path.size() >= 3 && path.compare(path.size() - 3, 3, ".gz") == 0;
    gzFile f = gzopen(path.c_str(), gz ? "wb6" : "wbT");
    if (f == nullptr) throw std::runtime_error("Cannot open " + path + " for writing: " + std::strerror(errno));
    file_ = f;
    buffer_.reserve(1 << 20);
}

OutputFile::~OutputFile() {
    if (file_ != nullptr) gzclose(static_cast<gzFile>(file_));
}

void OutputFile::write(const std::string& text) {
    buffer_ += text;
    if (buffer_.size() >= (1u << 20)) flush_buffer();
}

void OutputFile::flush_buffer() {
    if (buffer_.empty()) return;
    const int written = gzwrite(static_cast<gzFile>(file_), buffer_.data(), static_cast<unsigned>(buffer_.size()));
    if (written != static_cast<int>(buffer_.size())) {
        int errnum = 0;
        throw std::runtime_error("Error while writing " + path_ + ": " + gzerror(static_cast<gzFile>(file_), &errnum));
    }
    buffer_.clear();
}

void OutputFile::close() {
    if (file_ == nullptr) return;
    flush_buffer();
    const int status = gzclose(static_cast<gzFile>(file_));
    file_ = nullptr;
    if (status != Z_OK) throw std::runtime_error("Error while closing " + path_ + ".");
}

void Logger::info(const std::string& message) {
    if (quiet_) return;
    std::lock_guard<std::mutex> lock(mutex_);
    std::cerr << message << std::endl;
}

void Logger::warn(const std::string& message) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::cerr << message << std::endl;
}

}  // namespace charon
