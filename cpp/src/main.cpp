// charon: command line interface of the C++ port of Charon.jl.
//
//   charon [mcmc] [options] <dicefile> <outputfile>
//   charon [mcmc] [options] <basecountfile> <frequencyfile> <outputfile>
//   charon exactposterior [options] <dicefile> <outputfile>
//   charon exactposterior [options] <basecountfile> <frequencyfile> <outputfile>

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <exception>
#include <filesystem>
#include <iostream>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "charon/data.hpp"
#include "charon/io.hpp"
#include "charon/priors.hpp"
#include "charon/sampler.hpp"
#include "charon/version.hpp"

namespace {

using namespace charon;

const char* kUsage = R"(CHARON: Calculating HAplotype distribution to Resolve Original Number of individuals and demography.

Estimates the number of individuals n contributing to an environmental or forensic DNA
sample, the drift times tauC (anchor population) and tauA (sample population) and the
error rate eps, with the MCMC sampler of van Waaij et al. (2026).

Usage:
  charon [mcmc] [options] <dicefile> <outputfile>
  charon [mcmc] [options] <basecountfile> <frequencyfile> <outputfile>
  charon exactposterior [options] <dicefile> <outputfile>
  charon exactposterior [options] <basecountfile> <frequencyfile> <outputfile>

Input files may be gzipped. A DICE-2 file has four columns: ancestral reads, derived
reads, anchor frequency, count. A base count file has two columns: derived reads and
total reads; the frequency file has one column with the anchor frequency, one row per
row of the base count file. Columns may be separated by commas, tabs or spaces. The
output is CSV; it is gzipped when the output file name ends in .gz.

MCMC options (defaults as in runmcmc.jl):
  -c, --chains N          number of chains (default 4)
  -s, --steps N           number of samples per chain (default 100000)
      --n-prior SPEC      prior on n: uniform:A:B or weights:w1,w2,... (default uniform:1:10)
      --tauC-prior SPEC   prior on tauC: uniform:LO:HI, beta:A:B or beta:A:B:LO:HI
                          (default uniform:0:1)
      --tauA-prior SPEC   prior on tauA (default uniform:0:1)
      --eps-prior SPEC    prior on eps, support within [0, 0.5] (default uniform:0:0.5)
      --seed N            random seed (default: random); chains are reproducible for a
                          given seed, whatever the number of threads
      --messages N        progress message every N steps (default steps/100, 0: none)
      --scaling-messages  report changes of the adaptive proposal scale factors
      --raw-chains PREFIX also write every full chain, including the chains conditioned
                          on each n, to PREFIX_chain<k>.csv

exactposterior options (log-likelihood on a grid; defaults are the grid that the MCMC
sampler uses to find starting values):
      --n LIST            values of n, e.g. 1:10 or 1,2,4 (default 1:10)
      --tauC LIST         values of tauC: LO:HI:N (N equally spaced values) or a comma
                          separated list (default: 22 values 1/22, ..., 1)
      --tauA LIST         values of tauA (default as tauC)
      --eps LIST          values of eps in [0, 0.5) (default: 22 values 0.5/23, ..., 11/23)
      --messages N        progress message every N grid points (default: about 1%)

General options:
  -t, --threads N         number of threads (default: all cores)
      --header            the input files have a header line
      --no-header         the input files have no header line
                          (default: detect, a header is a first line that does not
                          start with a number)
  -q, --quiet             no progress messages
  -h, --help              show this help
  -V, --version           show the version

Example:
  charon --threads 4 IBS_ind4_basecount_IBS_ind4.gz IBS_ind4_freq_IBS.gz output.csv
  charon --threads 4 IBS_ind4.dice.gz output.csv
)";

struct Options {
    std::string command = "mcmc";
    std::vector<std::string> files;
    MCMCSettings mcmc;
    bool seed_given = false;
    bool messages_given = false;
    std::int64_t messages = 0;
    std::string raw_prefix;
    std::optional<bool> header;
    bool quiet = false;
    std::vector<int> grid_n;
    std::vector<double> grid_tauC, grid_tauA, grid_eps;
};

[[noreturn]] void usage_error(const std::string& message) {
    throw std::invalid_argument(message + "\nRun 'charon --help' for usage.");
}

bool is_plain_integer(const std::string& text, bool allow_sign) {
    if (text.empty()) return false;
    std::size_t i = 0;
    if (allow_sign && (text[0] == '-' || text[0] == '+')) i = 1;
    if (i >= text.size()) return false;
    for (; i < text.size(); ++i) {
        if (text[i] < '0' || text[i] > '9') return false;
    }
    return true;
}

long long parse_integer(const std::string& text, const std::string& option) {
    errno = 0;
    const long long v = is_plain_integer(text, true) ? std::strtoll(text.c_str(), nullptr, 10) : 0;
    if (!is_plain_integer(text, true) || errno == ERANGE) {
        usage_error("The value of " + option + " should be an integer, got '" + text + "'.");
    }
    return v;
}

std::uint64_t parse_seed(const std::string& text) {
    errno = 0;
    const unsigned long long v = is_plain_integer(text, false) ? std::strtoull(text.c_str(), nullptr, 10) : 0;
    if (!is_plain_integer(text, false) || errno == ERANGE) {
        usage_error("The seed should be an integer between 0 and 18446744073709551615, got '" + text + "'.");
    }
    return v;
}

/// Fail early, before hours of sampling, when the output file cannot be written, and
/// never overwrite an input file. Does not truncate an existing file.
void check_output_path(const std::string& path, const std::vector<std::string>& inputs) {
    std::error_code ec;
    for (const std::string& input : inputs) {
        if (std::filesystem::exists(path, ec) && std::filesystem::equivalent(path, input, ec)) {
            usage_error("The output file " + path + " is also an input file.");
        }
    }
    const bool existed = std::filesystem::exists(path, ec);
    std::FILE* f = std::fopen(path.c_str(), "ab");
    if (f == nullptr) {
        throw std::runtime_error("Cannot open " + path + " for writing: " + std::strerror(errno));
    }
    std::fclose(f);
    if (!existed) std::filesystem::remove(path, ec);
}

std::vector<int> parse_n_list(const std::string& text) {
    std::vector<int> ns;
    const std::size_t colon = text.find(':');
    if (colon != std::string::npos) {
        const long long a = parse_integer(text.substr(0, colon), "--n");
        const long long b = parse_integer(text.substr(colon + 1), "--n");
        if (a > b || b - a > 1000) usage_error("Invalid range for --n: '" + text + "'.");
        for (long long n = a; n <= b; ++n) ns.push_back(static_cast<int>(n));
    } else {
        for (double v : split_numbers(text, ',')) {
            if (v != std::floor(v) || std::fabs(v) > 1e6) usage_error("The values of --n should be integers.");
            ns.push_back(static_cast<int>(v));
        }
    }
    return ns;
}

std::vector<double> parse_grid(const std::string& text, const std::string& option) {
    if (text.find(':') != std::string::npos) {
        const std::vector<double> parts = split_numbers(text, ':');
        if (parts.size() != 3 || parts[2] < 1 || parts[2] != std::floor(parts[2]) || parts[2] > 1e7) {
            usage_error("The value of " + option + " should be LO:HI:N or a comma separated list, got '" + text + "'.");
        }
        const int count = static_cast<int>(parts[2]);
        std::vector<double> grid;
        for (int i = 0; i < count; ++i) {
            grid.push_back(count == 1 ? parts[0] : parts[0] + (parts[1] - parts[0]) * i / (count - 1.0));
        }
        return grid;
    }
    return split_numbers(text, ',');
}

Options parse_arguments(int argc, char** argv) {
    Options o;
    std::vector<std::string> args(argv + 1, argv + argc);
    std::size_t i = 0;
    if (!args.empty() && (args[0] == "mcmc" || args[0] == "exactposterior")) {
        o.command = args[0];
        i = 1;
    }
    auto value = [&](const std::string& option) -> std::string {
        if (i + 1 >= args.size()) usage_error("Option " + option + " needs a value.");
        return args[++i];
    };
    bool only_files = false;
    for (; i < args.size(); ++i) {
        std::string a = args[i];
        if (only_files || a.empty() || a[0] != '-' || a == "-") {
            o.files.push_back(a);
            continue;
        }
        // Allow --option=value for options that take a value.
        const std::size_t eq = a.find('=');
        if (a.rfind("--", 0) == 0 && eq != std::string::npos) {
            const std::string name = a.substr(0, eq);
            static const char* const flags[] = {"--help", "--version", "--scaling-messages", "--header",
                                                "--no-header", "--quiet"};
            for (const char* flag : flags) {
                if (name == flag) usage_error("Option " + name + " does not take a value.");
            }
            args.insert(args.begin() + static_cast<std::ptrdiff_t>(i) + 1, a.substr(eq + 1));
            a = name;
        }
        if (a == "--") {
            only_files = true;
        } else if (a == "-h" || a == "--help") {
            std::cout << kUsage;
            std::exit(0);
        } else if (a == "-V" || a == "--version") {
            std::cout << "charon " << kVersion << " (C++ port of Charon.jl " << kJuliaReferenceVersion << ")\n";
            std::exit(0);
        } else if (a == "-c" || a == "--chains") {
            const long long v = parse_integer(value(a), a);
            if (v < 1 || v > 100000) usage_error("The number of chains should be between 1 and 100000.");
            o.mcmc.nchains = static_cast<int>(v);
        } else if (a == "-s" || a == "--steps") {
            const long long v = parse_integer(value(a), a);
            if (v < 1 || v > 2000000000LL) usage_error("The number of steps should be between 1 and 2000000000.");
            o.mcmc.nsteps = static_cast<int>(v);
        } else if (a == "-t" || a == "--threads") {
            const long long v = parse_integer(value(a), a);
            if (v < 1 || v > 4096) usage_error("The number of threads should be between 1 and 4096.");
            o.mcmc.threads = static_cast<int>(v);
        } else if (a == "--n-prior") {
            o.mcmc.prioronn = DiscretePrior::parse(value(a));
        } else if (a == "--tauC-prior" || a == "--tauc-prior") {
            o.mcmc.prioronTauC = ContinuousPrior::parse(value(a));
        } else if (a == "--tauA-prior" || a == "--taua-prior") {
            o.mcmc.prioronTauA = ContinuousPrior::parse(value(a));
        } else if (a == "--eps-prior") {
            o.mcmc.prioronEps = ContinuousPrior::parse(value(a));
        } else if (a == "--seed") {
            o.mcmc.seed = parse_seed(value(a));
            o.seed_given = true;
        } else if (a == "--messages") {
            const long long v = parse_integer(value(a), a);
            if (v < 0) usage_error("The value of --messages should be non-negative.");
            o.messages = v;
            o.messages_given = true;
        } else if (a == "--scaling-messages") {
            o.mcmc.scalingmessages = true;
        } else if (a == "--raw-chains") {
            o.raw_prefix = value(a);
        } else if (a == "--header") {
            o.header = true;
        } else if (a == "--no-header") {
            o.header = false;
        } else if (a == "-q" || a == "--quiet") {
            o.quiet = true;
        } else if (a == "--n") {
            o.grid_n = parse_n_list(value(a));
        } else if (a == "--tauC" || a == "--tauc") {
            o.grid_tauC = parse_grid(value(a), a);
        } else if (a == "--tauA" || a == "--taua") {
            o.grid_tauA = parse_grid(value(a), a);
        } else if (a == "--eps") {
            o.grid_eps = parse_grid(value(a), a);
        } else {
            usage_error("Unknown option '" + a + "'.");
        }
    }
    if (o.files.size() != 2 && o.files.size() != 3) {
        usage_error(
            "You should provide either two or three file names: if you provide two, the first should be the DICE "
            "file and the second the output file. If you provide three, the first should be the base count file, "
            "the second the frequency file and the third the output file.");
    }
    return o;
}

std::string elapsed_since(std::chrono::steady_clock::time_point start) {
    const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.1f s", seconds);
    return buf;
}

PreparedData load_data(const Options& o, Logger& log) {
    ReadSummary summary;
    Observations obs;
    if (o.files.size() == 2) {
        log.info("Reading the DICE file " + o.files[0] + ".");
        obs = read_dice_file(o.files[0], o.header, &summary);
    } else {
        log.info("Reading the base count file " + o.files[0] + " and the frequency file " + o.files[1] + ".");
        obs = read_basecount_and_frequency_files(o.files[0], o.files[1], o.header, o.header, &summary);
    }
    PreparedData data = prepare_data(obs);
    log.info("Read " + std::to_string(summary.rows) + " rows (" + std::to_string(summary.loci) + " loci); " +
             std::to_string(static_cast<long long>(data.total_count())) +
             " loci have positive coverage and an anchor frequency strictly between 0 and 1, in " +
             std::to_string(data.num_loci()) + " distinct (coverage, derived reads, frequency) combinations.");
    return data;
}

int run(const Options& o) {
    Logger log(o.quiet);
    const auto start = std::chrono::steady_clock::now();
    const std::string& outputfile = o.files.back();
    const std::vector<std::string> inputs(o.files.begin(), o.files.end() - 1);
    check_output_path(outputfile, inputs);
    if (!o.raw_prefix.empty() && o.command != "exactposterior") {
        for (int c = 1; c <= o.mcmc.nchains; ++c) {
            check_output_path(o.raw_prefix + "_chain" + std::to_string(c) + ".csv", inputs);
        }
    }
    PreparedData data = load_data(o, log);

    if (o.command == "exactposterior") {
        const std::vector<int> ns = o.grid_n.empty() ? parse_n_list("1:10") : o.grid_n;
        std::vector<double> defaulttau, defaulteps;
        for (int i = 1; i <= 22; ++i) {
            defaulttau.push_back(i / 22.0);
            defaulteps.push_back(0.5 * i / 23.0);
        }
        const auto& tauC = o.grid_tauC.empty() ? defaulttau : o.grid_tauC;
        const auto& tauA = o.grid_tauA.empty() ? defaulttau : o.grid_tauA;
        const auto& eps = o.grid_eps.empty() ? defaulteps : o.grid_eps;
        const std::size_t total = ns.size() * tauC.size() * tauA.size() * eps.size();
        const std::int64_t messages =
            o.quiet ? 0 : (o.messages_given ? o.messages : static_cast<std::int64_t>(total / 100));
        log.info("Calculating the log-likelihood on " + std::to_string(total) + " grid points.");
        const GridResult grid = exactposterior(data, ns, tauC, tauA, eps, o.mcmc.threads, messages, log);
        OutputFile out(outputfile);
        write_grid(grid, out);
        out.close();
        log.info("I am finished after " + elapsed_since(start) + ". I saved the results in " + outputfile + ".");
        return 0;
    }

    MCMCSettings settings = o.mcmc;
    if (!o.seed_given) settings.seed = (static_cast<std::uint64_t>(std::random_device{}()) << 32) ^ std::random_device{}();
    settings.messages = o.quiet ? 0 : (o.messages_given ? o.messages : -1);
    log.info("Priors: n ~ " + settings.prioronn.describe() + ", τC ~ " + settings.prioronTauC.describe() +
             ", τA ~ " + settings.prioronTauA.describe() + ", ϵ ~ " + settings.prioronEps.describe() + ".");
    log.info("Seed: " + std::to_string(settings.seed) + ".");
    const std::vector<Chain> chains = mcmc_sampler(data, settings, log);

    log.info("I am unpacking the posterior and writing it to " + outputfile + ".");
    OutputFile out(outputfile);
    write_unpacked_posterior(chains, out);
    out.close();
    if (!o.raw_prefix.empty()) {
        for (std::size_t c = 0; c < chains.size(); ++c) {
            const std::string path = o.raw_prefix + "_chain" + std::to_string(c + 1) + ".csv";
            OutputFile raw(path);
            write_raw_chain(chains[c], raw);
            raw.close();
        }
        log.info("I saved the full chains as " + o.raw_prefix + "_chain<k>.csv.");
    }
    log.info("I am finished after " + elapsed_since(start) + ". I saved the results in " + outputfile + ".");
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        return run(parse_arguments(argc, argv));
    } catch (const std::invalid_argument& e) {
        std::cerr << "charon: error: " << e.what() << std::endl;
        return 2;
    } catch (const std::exception& e) {
        std::cerr << "charon: error: " << e.what() << std::endl;
        return 1;
    }
}
