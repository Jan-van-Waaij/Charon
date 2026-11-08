#include "charon/data.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>

#include <zlib.h>

namespace charon {

namespace {

/// Buffered line reader for plain or gzipped files (zlib reads plain files transparently).
class LineReader {
public:
    explicit LineReader(const std::string& path) : path_(path), buffer_(1 << 22) {
        file_ = gzopen(path.c_str(), "rb");
        if (file_ == nullptr) throw std::runtime_error("Cannot open file " + path + ": " + std::strerror(errno));
        gzbuffer(file_, 1 << 20);
    }
    ~LineReader() {
        if (file_ != nullptr) gzclose(file_);
    }
    LineReader(const LineReader&) = delete;
    LineReader& operator=(const LineReader&) = delete;

    /// Next line without its line ending; the view stays valid until the next call.
    bool getline(std::string_view& line) {
        for (;;) {
            const char* start = buffer_.data() + begin_;
            const void* newline = std::memchr(start, '\n', end_ - begin_);
            if (newline != nullptr) {
                const std::size_t length = static_cast<const char*>(newline) - start;
                line = trim_cr(std::string_view(start, length));
                begin_ += length + 1;
                ++line_number_;
                return true;
            }
            if (eof_) {
                if (begin_ == end_) return false;
                line = trim_cr(std::string_view(start, end_ - begin_));
                begin_ = end_;
                ++line_number_;
                return true;
            }
            fill();
        }
    }

    std::uint64_t line_number() const { return line_number_; }
    const std::string& path() const { return path_; }

private:
    static std::string_view trim_cr(std::string_view s) {
        if (!s.empty() && s.back() == '\r') s.remove_suffix(1);
        return s;
    }

    void fill() {
        // Move the incomplete line to the front, grow the buffer for very long lines.
        const std::size_t remaining = end_ - begin_;
        if (begin_ > 0 && remaining > 0) std::memmove(buffer_.data(), buffer_.data() + begin_, remaining);
        begin_ = 0;
        end_ = remaining;
        if (end_ == buffer_.size()) buffer_.resize(buffer_.size() * 2);
        const unsigned want = static_cast<unsigned>(std::min<std::size_t>(buffer_.size() - end_, 1u << 30));
        const int got = gzread(file_, buffer_.data() + end_, want);
        if (got < 0) {
            int errnum = 0;
            const char* message = gzerror(file_, &errnum);
            throw std::runtime_error("Error while reading " + path_ + ": " + message);
        }
        // A truncated gzip file sets Z_BUF_ERROR on the read that returns the last data.
        // zlib >= 1.3.2 clears that error on the next call, so check after every read.
        int errnum = 0;
        const char* message = gzerror(file_, &errnum);
        if (errnum != Z_OK && errnum != Z_STREAM_END) {
            throw std::runtime_error("Error while reading " + path_ + ": " + message +
                                     " (is the gzip file truncated?)");
        }
        end_ += static_cast<std::size_t>(got);
        if (got == 0) eof_ = true;
    }

    std::string path_;
    gzFile file_ = nullptr;
    std::vector<char> buffer_;
    std::size_t begin_ = 0;
    std::size_t end_ = 0;
    bool eof_ = false;
    std::uint64_t line_number_ = 0;
};

inline bool is_separator(char c) {
    return c == ' ' || c == '\t' || c == ',' || c == ';' || c == '|';
}

/// Split a line into at most max_fields fields; returns the number of fields found
/// (which can exceed max_fields, in which case the extra fields are not stored).
std::size_t split_fields(std::string_view line, std::string_view* fields, std::size_t max_fields) {
    std::size_t count = 0;
    std::size_t i = 0;
    const std::size_t n = line.size();
    while (i < n) {
        while (i < n && is_separator(line[i])) ++i;
        if (i >= n) break;
        const std::size_t start = i;
        while (i < n && !is_separator(line[i])) ++i;
        if (count < max_fields) fields[count] = line.substr(start, i - start);
        ++count;
    }
    return count;
}

bool is_blank(std::string_view line) {
    for (char c : line) {
        if (c != ' ' && c != '\t') return false;
    }
    return true;
}

bool parse_double(std::string_view token, double& value) {
    char buf[128];
    if (token.empty() || token.size() >= sizeof(buf)) return false;
    std::memcpy(buf, token.data(), token.size());
    buf[token.size()] = '\0';
    char* end = nullptr;
    errno = 0;
    value = std::strtod(buf, &end);
    return end == buf + token.size() && errno != ERANGE;
}

bool parse_int(std::string_view token, std::int64_t& value) {
    char buf[64];
    if (token.empty() || token.size() >= sizeof(buf)) return false;
    std::memcpy(buf, token.data(), token.size());
    buf[token.size()] = '\0';
    char* end = nullptr;
    errno = 0;
    const long long v = std::strtoll(buf, &end, 10);
    if (end == buf + token.size() && errno == 0) {
        value = v;
        return true;
    }
    // Accept integral floating point notation such as "3.0".
    double d = 0.0;
    if (parse_double(token, d) && std::isfinite(d) && d == std::floor(d) && std::fabs(d) < 9.0e15) {
        value = static_cast<std::int64_t>(d);
        return true;
    }
    return false;
}

[[noreturn]] void fail(const LineReader& reader, const std::string& message) {
    throw std::runtime_error(reader.path() + ", line " + std::to_string(reader.line_number()) + ": " + message);
}

/// The next non-empty line split into exactly `expected` fields. Handles the header.
class Table {
public:
    Table(const std::string& path, std::size_t expected, std::optional<bool> header, const char* description)
        : reader_(path), expected_(expected), description_(description) {
        std::string_view line;
        if (!next_nonblank(line)) {
            throw std::runtime_error("The " + std::string(description_) + " " + path + " is empty.");
        }
        if (line.size() >= 3 && line.compare(0, 3, "\xEF\xBB\xBF") == 0) line.remove_prefix(3);  // UTF-8 BOM
        std::string_view first[8];
        const std::size_t nfields = split_fields(line, first, 8);
        double dummy = 0.0;
        has_header_ = header.has_value() ? *header : !(nfields > 0 && parse_double(first[0], dummy));
        if (!has_header_) {
            pending_ = line;
            has_pending_ = true;
        }
    }

    bool has_header() const { return has_header_; }

    bool next(std::string_view* fields) {
        std::string_view line;
        if (has_pending_) {
            line = pending_;
            has_pending_ = false;
        } else if (!next_nonblank(line)) {
            return false;
        }
        const std::size_t nfields = split_fields(line, fields, expected_);
        if (nfields != expected_) {
            fail(reader_, "expected " + std::to_string(expected_) + " columns in the " + description_ + ", found " +
                              std::to_string(nfields) + ".");
        }
        return true;
    }

    std::int64_t integer(std::string_view token, const char* column) const {
        std::int64_t v = 0;
        if (!parse_int(token, v)) {
            fail(reader_, "could not read the " + std::string(column) + " '" + std::string(token) + "' as an integer.");
        }
        return v;
    }

    double real(std::string_view token, const char* column) const {
        double v = 0.0;
        if (!parse_double(token, v)) {
            fail(reader_, "could not read the " + std::string(column) + " '" + std::string(token) + "' as a number.");
        }
        return v;
    }

    const LineReader& reader() const { return reader_; }

private:
    bool next_nonblank(std::string_view& line) {
        while (reader_.getline(line)) {
            if (!is_blank(line)) return true;
        }
        return false;
    }

    LineReader reader_;
    std::size_t expected_;
    const char* description_;
    bool has_header_ = false;
    bool has_pending_ = false;
    std::string_view pending_;  // points into the reader's buffer, consumed before the next read
};

void check_row(const LineReader& reader, std::int64_t coverage, std::int64_t derived, double frequency,
               std::int64_t count) {
    if (derived < 0) fail(reader, "the number of derived reads should be non-negative.");
    if (derived > coverage) fail(reader, "the number of derived reads should be at most the coverage.");
    if (!(frequency >= 0.0 && frequency <= 1.0)) fail(reader, "the frequency should be between 0.0 and 1.0.");
    if (count < 0) fail(reader, "the count should be non-negative.");
}

struct LocusKey {
    std::int64_t coverage;
    std::int64_t derived;
    std::uint64_t frequency_bits;
    bool operator==(const LocusKey& o) const {
        return coverage == o.coverage && derived == o.derived && frequency_bits == o.frequency_bits;
    }
};

struct LocusKeyHash {
    std::size_t operator()(const LocusKey& k) const {
        std::uint64_t h = k.frequency_bits * 0x9E3779B97F4A7C15ULL;
        h ^= static_cast<std::uint64_t>(k.coverage) * 0xC2B2AE3D27D4EB4FULL + (h << 6) + (h >> 2);
        h ^= static_cast<std::uint64_t>(k.derived) * 0x165667B19E3779F9ULL + (h << 6) + (h >> 2);
        return static_cast<std::size_t>(h ^ (h >> 29));
    }
};

std::uint64_t bits_of(double x) {
    std::uint64_t b;
    std::memcpy(&b, &x, sizeof b);
    return b;
}

}  // namespace

void Observations::push_back(std::int64_t coverage, std::int64_t derived, double frequency, std::int64_t count) {
    coverages.push_back(coverage);
    derivedreads.push_back(derived);
    frequencies.push_back(frequency);
    counts.push_back(count);
}

Observations read_dice_file(const std::string& path, std::optional<bool> header, ReadSummary* summary) {
    Table table(path, 4, header, "DICE file");
    Observations obs;
    std::string_view f[4];
    std::uint64_t rows = 0;
    std::uint64_t loci = 0;
    while (table.next(f)) {
        const std::int64_t ancestral = table.integer(f[0], "number of ancestral reads");
        const std::int64_t derived = table.integer(f[1], "number of derived reads");
        const double frequency = table.real(f[2], "frequency");
        const std::int64_t count = table.integer(f[3], "count");
        if (ancestral < 0) fail(table.reader(), "the number of ancestral reads should be non-negative.");
        if (derived > kMaxCoverage || ancestral > kMaxCoverage - std::max<std::int64_t>(derived, 0)) {
            fail(table.reader(), "the coverage is larger than " + std::to_string(kMaxCoverage) + ".");
        }
        check_row(table.reader(), ancestral + derived, derived, frequency, count);
        obs.push_back(ancestral + derived, derived, frequency, count);
        ++rows;
        loci += static_cast<std::uint64_t>(count);
    }
    if (rows == 0) throw std::runtime_error("The DICE file " + path + " contains no data.");
    if (summary != nullptr) {
        summary->rows = rows;
        summary->loci = loci;
        summary->header_first = table.has_header();
    }
    return obs;
}

Observations read_basecount_and_frequency_files(const std::string& basecountpath, const std::string& frequencypath,
                                                std::optional<bool> header_basecount,
                                                std::optional<bool> header_frequency, ReadSummary* summary) {
    Table basecounts(basecountpath, 2, header_basecount, "base count file");
    Table frequencies(frequencypath, 1, header_frequency, "frequency file");
    std::unordered_map<LocusKey, std::int64_t, LocusKeyHash> aggregated;
    std::vector<LocusKey> order;  // first-seen order, keeps the output deterministic
    std::string_view b[2];
    std::string_view f[1];
    std::uint64_t rows = 0;
    for (;;) {
        const bool more_b = basecounts.next(b);
        const bool more_f = frequencies.next(f);
        if (!more_b && !more_f) break;
        if (more_b != more_f) {
            throw std::runtime_error("The base count file " + basecountpath + " and the frequency file " +
                                     frequencypath + " should have the same number of rows, but the " +
                                     (more_b ? "frequency" : "base count") + " file ends after " +
                                     std::to_string(rows) + " data rows.");
        }
        const std::int64_t derived = basecounts.integer(b[0], "number of derived reads");
        const std::int64_t coverage = basecounts.integer(b[1], "total number of reads");
        const double frequency = frequencies.real(f[0], "frequency");
        ++rows;
        // As in Charon.jl, rows with zero coverage or with an anchor frequency that is not
        // strictly between 0 and 1 (including NaN) are removed before the data are checked.
        if (!(coverage > 0 && frequency > 0.0 && frequency < 1.0)) continue;
        check_row(basecounts.reader(), coverage, derived, 0.5, 1);
        if (coverage > kMaxCoverage) {
            fail(basecounts.reader(), "the coverage is larger than " + std::to_string(kMaxCoverage) + ".");
        }
        const LocusKey key{coverage, derived, bits_of(frequency)};
        auto [it, inserted] = aggregated.try_emplace(key, 0);
        if (inserted) order.push_back(key);
        ++it->second;
    }
    if (rows == 0) throw std::runtime_error("The base count file " + basecountpath + " contains no data.");
    if (order.empty()) {
        throw std::runtime_error(
            "At least one row should have positive coverage and an anchor frequency strictly between 0 and 1.");
    }
    Observations obs;
    for (const LocusKey& key : order) {
        double frequency;
        std::memcpy(&frequency, &key.frequency_bits, sizeof frequency);
        obs.push_back(key.coverage, key.derived, frequency, aggregated[key]);
    }
    if (summary != nullptr) {
        summary->rows = rows;
        summary->loci = rows;
        summary->header_first = basecounts.has_header();
        summary->header_second = frequencies.has_header();
    }
    return obs;
}

PreparedData prepare_data(const Observations& obs) {
    const std::size_t n = obs.size();
    if (obs.derivedreads.size() != n || obs.frequencies.size() != n || obs.counts.size() != n) {
        throw std::invalid_argument("coverages, derivedreads, frequencies and counts should be of the same length.");
    }
    if (n == 0) throw std::invalid_argument("The data is empty.");
    // The checks of Charon.jl's MCMCsampler.
    if (*std::max_element(obs.coverages.begin(), obs.coverages.end()) <= 0) {
        throw std::invalid_argument("At least some coverages should be positive.");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (obs.derivedreads[i] < 0) throw std::invalid_argument("All reads should be non-negative.");
        if (obs.derivedreads[i] > obs.coverages[i]) {
            throw std::invalid_argument("All reads should be at most the coverage.");
        }
        if (!(obs.frequencies[i] >= 0.0 && obs.frequencies[i] <= 1.0)) {
            throw std::invalid_argument("The frequencies should be between 0.0 and 1.0.");
        }
        if (obs.counts[i] < 0) throw std::invalid_argument("All counts should be non-negative.");
        if (obs.coverages[i] > kMaxCoverage) {
            throw std::invalid_argument("All coverages should be at most " + std::to_string(kMaxCoverage) + ".");
        }
    }

    // Keep loci with positive coverage and anchor frequency strictly between 0 and 1.
    std::vector<std::size_t> keep;
    keep.reserve(n);
    for (std::size_t i = 0; i < n; ++i) {
        if (obs.coverages[i] > 0 && obs.frequencies[i] > 0.0 && obs.frequencies[i] < 1.0 && obs.counts[i] > 0) {
            keep.push_back(i);
        }
    }
    if (keep.empty()) {
        throw std::invalid_argument(
            "At least one locus with positive coverage should have an anchor frequency strictly between 0 and 1.");
    }
    std::sort(keep.begin(), keep.end(), [&](std::size_t a, std::size_t b) {
        if (obs.coverages[a] != obs.coverages[b]) return obs.coverages[a] < obs.coverages[b];
        if (obs.derivedreads[a] != obs.derivedreads[b]) return obs.derivedreads[a] < obs.derivedreads[b];
        return obs.frequencies[a] < obs.frequencies[b];
    });

    PreparedData data;
    for (std::size_t idx : keep) {
        const std::int64_t R = obs.coverages[idx];
        const std::int64_t d = obs.derivedreads[idx];
        const double y = obs.frequencies[idx];
        const bool same_pair = !data.pair_coverage.empty() && data.pair_coverage.back() == R &&
                               data.pair_derived.back() == d;
        if (same_pair && data.frequencies.back() == y) {
            data.counts.back() += static_cast<double>(obs.counts[idx]);  // identical locus: aggregate
            continue;
        }
        if (!same_pair) {
            data.pair_coverage.push_back(R);
            data.pair_derived.push_back(d);
        }
        data.pair.push_back(static_cast<std::int32_t>(data.pair_coverage.size() - 1));
        data.frequencies.push_back(y);
        data.counts.push_back(static_cast<double>(obs.counts[idx]));
    }
    return data;
}

}  // namespace charon
