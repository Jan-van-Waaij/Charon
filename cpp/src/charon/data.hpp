#pragma once

// Reading and preparing the input data (the data handling of Charon.jl).
//
// Two input formats are supported, plain or gzipped:
//  * DICE-2: four columns, ancestral reads, derived reads, anchor frequency, count.
//  * a base count file with two columns (derived reads, total reads) together with a
//    frequency file with one column (anchor frequency), with the same number of rows.
// Columns may be separated by commas, tabs, spaces, semicolons or '|'.

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "charon/model.hpp"

namespace charon {

/// Largest coverage accepted (a sanity limit; real data have far lower coverage).
inline constexpr std::int64_t kMaxCoverage = 10000000;

/// Raw observations: counts[i] loci with coverage coverages[i], derivedreads[i]
/// derived reads and anchor frequency frequencies[i].
struct Observations {
    std::vector<std::int64_t> coverages;
    std::vector<std::int64_t> derivedreads;
    std::vector<double> frequencies;
    std::vector<std::int64_t> counts;

    std::size_t size() const { return coverages.size(); }
    void push_back(std::int64_t coverage, std::int64_t derived, double frequency, std::int64_t count);
};

/// Summary of what was read, for reporting.
struct ReadSummary {
    std::uint64_t rows = 0;          ///< data rows read (excluding a header)
    std::uint64_t loci = 0;          ///< number of loci represented (sum of counts)
    bool header_first = false;       ///< first (or only) file had a header
    bool header_second = false;      ///< frequency file had a header
};

/// Read a DICE-2 file (plain or gzipped). header: true/false, or nullopt to auto-detect
/// (the file has a header when the first field of the first line is not a number).
Observations read_dice_file(const std::string& path, std::optional<bool> header, ReadSummary* summary = nullptr);

/// Read a base count file (derived reads, total reads) and a frequency file (one
/// column) with the same number of rows. As in Charon.jl, rows with zero coverage or
/// with an anchor frequency outside (0, 1) are dropped before the remaining rows are
/// checked, and rows with identical data are aggregated.
Observations read_basecount_and_frequency_files(const std::string& basecountpath, const std::string& frequencypath,
                                                std::optional<bool> header_basecount,
                                                std::optional<bool> header_frequency,
                                                ReadSummary* summary = nullptr);

/// Check the data as Charon.jl does, remove loci with zero coverage or with anchor
/// frequency 0 or 1 (and loci with count 0), and aggregate identical loci.
/// Throws std::invalid_argument when the data cannot be used.
PreparedData prepare_data(const Observations& obs);

}  // namespace charon
