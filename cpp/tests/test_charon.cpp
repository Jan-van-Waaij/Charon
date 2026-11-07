// Unit tests for the C++ port of Charon.jl. Several tests are ports of the tests in
// Charon.jl's test/runtests.jl; the reference log-likelihoods were computed with
// Charon.jl 0.3.2 (Julia 1.10.5) on tests/data/small.dice.csv.

#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <zlib.h>

#include "charon/data.hpp"
#include "charon/io.hpp"
#include "charon/model.hpp"
#include "charon/priors.hpp"
#include "charon/rng.hpp"
#include "charon/sampler.hpp"

using namespace charon;
namespace fs = std::filesystem;

namespace {

int g_checks = 0;
int g_failures = 0;

#define CHECK(cond)                                                                       \
    do {                                                                                  \
        ++g_checks;                                                                       \
        if (!(cond)) {                                                                    \
            ++g_failures;                                                                 \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n";    \
        }                                                                                 \
    } while (0)

#define CHECK_CLOSE(a, b, rtol)                                                                              \
    do {                                                                                                     \
        ++g_checks;                                                                                          \
        const double va_ = (a), vb_ = (b);                                                                   \
        if (!(std::fabs(va_ - vb_) <= (rtol) * std::max(std::fabs(va_), std::fabs(vb_)) || va_ == vb_)) {    \
            ++g_failures;                                                                                    \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK_CLOSE failed: " #a " = " << va_ << ", " #b " = " \
                      << vb_ << "\n";                                                                        \
        }                                                                                                    \
    } while (0)

#define CHECK_THROWS(expr, Exception)                                                                 \
    do {                                                                                              \
        ++g_checks;                                                                                   \
        bool thrown_ = false;                                                                         \
        try {                                                                                         \
            (void)(expr);                                                                             \
        } catch (const Exception&) {                                                                  \
            thrown_ = true;                                                                           \
        }                                                                                             \
        if (!thrown_) {                                                                               \
            ++g_failures;                                                                             \
            std::cerr << __FILE__ << ":" << __LINE__ << ": expected " #Exception " from " #expr "\n"; \
        }                                                                                             \
    } while (0)

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

std::string data_file(const std::string& name) { return std::string(CHARON_TEST_DATA_DIR) + "/" + name; }

fs::path temp_dir() {
    static const fs::path dir = [] {
        fs::path d = fs::temp_directory_path() /
                     ("charon_tests_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
        fs::create_directories(d);
        return d;
    }();
    return dir;
}

void write_text(const fs::path& path, const std::string& text) {
    if (path.extension() == ".gz") {
        gzFile f = gzopen(path.string().c_str(), "wb");
        gzwrite(f, text.data(), static_cast<unsigned>(text.size()));
        gzclose(f);
    } else {
        std::ofstream(path, std::ios::binary) << text;
    }
}

std::vector<std::vector<std::string>> read_csv(const fs::path& path) {
    std::vector<std::vector<std::string>> rows;
    gzFile f = gzopen(path.string().c_str(), "rb");
    std::string content;
    char buf[1 << 16];
    int got;
    while ((got = gzread(f, buf, sizeof buf)) > 0) content.append(buf, static_cast<std::size_t>(got));
    gzclose(f);
    std::istringstream in(content);
    std::string line;
    while (std::getline(in, line)) {
        std::vector<std::string> fields;
        std::istringstream ls(line);
        std::string field;
        while (std::getline(ls, field, ',')) fields.push_back(field);
        rows.push_back(fields);
    }
    return rows;
}

using LMat = Eigen::Matrix<long double, Eigen::Dynamic, Eigen::Dynamic>;

/// Reference matrix exponential: scaling and squaring with a Taylor series in long double.
LMat expm_reference(const Mat& A, double t) {
    LMat X = A.cast<long double>() * static_cast<long double>(t);
    int squarings = 0;
    long double norm = X.cwiseAbs().colwise().sum().maxCoeff();
    while (norm > 0.25L) {
        X /= 2;
        norm /= 2;
        ++squarings;
    }
    LMat result = LMat::Identity(X.rows(), X.cols());
    LMat term = result;
    for (int k = 1; k <= 30; ++k) {
        term = term * X / static_cast<long double>(k);
        result += term;
    }
    for (int i = 0; i < squarings; ++i) result = result * result;
    return result;
}

long double binomial_pdf_reference(int R, long double p, int d) {
    long double c = 1;
    for (int i = 1; i <= d; ++i) c = c * (R - d + i) / i;
    return c * std::pow(p, static_cast<long double>(d)) * std::pow(1 - p, static_cast<long double>(R - d));
}

PreparedData small_data() { return prepare_data(read_dice_file(data_file("small.dice.csv"), std::nullopt)); }

// ------------------------------------------------------------------------- tests

void test_makeQ() {
    // Ports of the makeQ and makeQꜜ tests.
    Mat expectedQ1(3, 3);
    expectedQ1 << 0, 1, 0, 0, -1, 0, 0, 1, 0;
    CHECK(makeQ(1) == expectedQ1);
    Mat expectedQdown1(3, 3);
    expectedQdown1 << 0, 3, 0, 0, -2, 1, 0, 1, -2;
    CHECK(makeQdown(1) == expectedQdown1);
    for (int n = 1; n <= 10; ++n) {
        const Mat Q = makeQ(n);
        const Mat Qd = makeQdown(n);
        for (int k = 0; k <= 2 * n; ++k) {
            for (int L = 0; L <= 2 * n; ++L) {
                double q = 0.0, qd = 0.0;
                if (L == k - 1) {
                    q = qd = 0.5 * k * (k - 1);
                } else if (L == k) {
                    q = -k * (2.0 * n - k);
                    qd = -k * (2.0 * n - k + 1);
                } else if (L == k + 1) {
                    q = 0.5 * (2.0 * n - k) * (2.0 * n - k - 1);
                    qd = 0.5 * (2.0 * n - k + 1) * (2.0 * n - k);
                }
                CHECK(Q(k, L) == q);
                CHECK(Qd(k, L) == qd);
            }
        }
    }
}

void test_eigen_expansion() {
    for (int n = 1; n <= 10; ++n) {
        for (const Mat& A : {makeQ(n), makeQdown(n)}) {
            const EigenExpansion F(A);
            const Mat back = F.P * F.D.asDiagonal() * F.Pinv;
            CHECK((back - A).cwiseAbs().maxCoeff() <= 1e-9 * A.cwiseAbs().maxCoeff());
            CHECK((F.exp(0.0) - Mat::Identity(A.rows(), A.cols())).cwiseAbs().maxCoeff() < 1e-8);
            for (double t : {0.01, 0.25, 1.0}) {
                const Mat E = F.exp(t);
                const LMat R = expm_reference(A, t);
                const double err = static_cast<double>((E.cast<long double>() - R).cwiseAbs().maxCoeff());
                // The eigenvectors of Q and Qdown are ill-conditioned for larger n
                // (cond(P) ~ 1e6 for n = 10), which limits the accuracy to ~1e-9.
                CHECK(err <= 1e-8 * static_cast<double>(R.cwiseAbs().maxCoeff()));
            }
        }
        // calcmatrix: e^{tauA Q} e^{tauC Qdown}
        const EigenExpansion Q(makeQ(n)), Qd(makeQdown(n));
        CHECK((calcmatrix(0.0, 0.0, Qd, Q) - Mat::Identity(2 * n + 1, 2 * n + 1)).cwiseAbs().maxCoeff() < 1e-8);
        const Mat M = calcmatrix(0.3, 0.2, Qd, Q);
        CHECK((M - Q.exp(0.2) * Qd.exp(0.3)).cwiseAbs().maxCoeff() == 0.0);
    }
    // A general symmetric matrix, and a matrix with complex eigenvalues.
    Mat S(3, 3);
    S << 2, 1, 0, 1, 3, 1, 0, 1, 4;
    const EigenExpansion FS(S);
    const LMat RS = expm_reference(S, 0.5);
    CHECK(static_cast<double>((FS.exp(0.5).cast<long double>() - RS).cwiseAbs().maxCoeff()) < 1e-12);
    // A badly scaled matrix, which the balancing step rescales.
    Eigen::Vector3d scaling(1.0, 1e4, 1e-4);
    const Mat badly = scaling.cwiseInverse().asDiagonal() * S * scaling.asDiagonal();
    const EigenExpansion FB(badly);
    const LMat RB = expm_reference(badly, 0.5);
    CHECK(static_cast<double>((FB.exp(0.5).cast<long double>() - RB).cwiseAbs().maxCoeff()) <
          1e-12 * static_cast<double>(RB.cwiseAbs().maxCoeff()));
    Mat rotation(2, 2);
    rotation << 0, -1, 1, 0;  // eigenvalues ±i
    CHECK_THROWS(EigenExpansion(rotation), std::runtime_error);
}

void test_binomial() {
    const std::vector<double> c = binomial_coefficients(20);
    CHECK(c[0] == 1 && c[1] == 20 && c[10] == 184756 && c[20] == 1);
    CHECK(binomial_coefficients(2) == std::vector<double>({1, 2, 1}));
    for (int R : {1, 2, 7, 30, 100}) {
        for (double p : {0.001, 0.1, 0.5, 0.73, 0.999}) {
            for (int d = 0; d <= R; ++d) {
                const long double ref = binomial_pdf_reference(R, p, d);
                if (ref < 1e-250L) continue;
                CHECK_CLOSE(binomial_pdf(R, p, d), static_cast<double>(ref), 1e-12);
            }
        }
    }
    CHECK(binomial_pdf(5, 0.0, 0) == 1.0);
    CHECK(binomial_pdf(5, 0.0, 1) == 0.0);
    CHECK(binomial_pdf(5, 1.0, 5) == 1.0);
    CHECK(binomial_pdf(5, 1.0, 4) == 0.0);
    CHECK(binomial_pdf(5, 0.3, 6) == 0.0);
    CHECK(binomial_pdf(5, 0.3, -1) == 0.0);
}

void test_loglik_reference_values() {
    // Values from Charon.exactposterior (Charon.jl 0.3.2) on tests/data/small.dice.csv.
    struct Ref {
        int n;
        double tauC, tauA, eps, loglik, rtol;
    };
    const Ref refs[] = {
        {1, 0.0, 0.0, 0.0, -8163.3520788234755, 1e-13},
        {1, 0.05, 0.2, 0.01, -6564.313085557253, 1e-13},
        {2, 0.3, 0.7, 0.01, -6469.3587446804395, 1e-13},
        {3, 0.05, 0.2, 0.01, -5854.51423556733, 1e-12},
        {3, 0.9, 0.0, 0.3, -12927.585079424049, 1e-12},
        {5, 0.3, 0.2, 0.01, -6010.457620026086, 1e-10},
    };
    const PreparedData data = small_data();
    for (const Ref& r : refs) {
        const ConditionalModel model(r.n, data);
        auto ws = model.make_workspace();
        CHECK_CLOSE(model.loglik(r.tauC, r.tauA, r.eps, ws), r.loglik, r.rtol);
    }
}

void test_loglik_without_drift() {
    // Port of the logprobderivedreads! test: with tauC = tauA = 0 the likelihood of a locus
    // is sum_k Binomial(R, p_k)(d) Binomial(2n, y)(k).
    const PreparedData data = small_data();
    for (int n = 1; n <= 10; ++n) {
        const ConditionalModel model(n, data);
        auto ws = model.make_workspace();
        CHECK(model.loglik(-0.1, 0.2, 0.1, ws) == kNegInf);
        CHECK(model.loglik(0.3, -0.1, 0.1, ws) == kNegInf);
        CHECK(model.loglik(0.3, 0.2, -0.1, ws) == kNegInf);
        CHECK(model.loglik(0.3, 0.2, 1.1, ws) == kNegInf);
        for (double eps : {0.1, 0.2, 0.3, 0.4}) {
            long double expected = 0;
            for (std::size_t i = 0; i < data.num_loci(); ++i) {
                const int R = static_cast<int>(data.pair_coverage[data.pair[i]]);
                const int d = static_cast<int>(data.pair_derived[data.pair[i]]);
                long double prob = 0;
                for (int k = 0; k <= 2 * n; ++k) {
                    const long double pk = k * (1.0L - eps) / (2 * n) + (2 * n - k) * static_cast<long double>(eps) / (2 * n);
                    prob += binomial_pdf_reference(R, pk, d) * binomial_pdf_reference(2 * n, data.frequencies[i], k);
                }
                expected += data.counts[i] * std::log(prob);
            }
            // n >= 8 is affected by the conditioning of the eigenvectors (see the README).
            const double rtol = n <= 5 ? 1e-9 : 1e-5;
            CHECK_CLOSE(model.loglik(0.0, 0.0, eps, ws), static_cast<double>(expected), rtol);
        }
    }
}

void test_prepare_data() {
    // Port of the filtervectorsandapplycountmap test.
    Observations obs;
    obs.push_back(3, 1, 0.25, 1);
    obs.push_back(3, 1, 0.25, 1);
    obs.push_back(2, 2, 0.5, 1);
    obs.push_back(3, 3, 0.75, 1);
    obs.push_back(4, 4, 0.875, 1);
    obs.push_back(5, 1, 0.0, 7);   // frequency 0: removed
    obs.push_back(5, 1, 1.0, 7);   // frequency 1: removed
    obs.push_back(0, 0, 0.5, 7);   // coverage 0: removed
    obs.push_back(6, 2, 0.3, 0);   // count 0: removed
    const PreparedData data = prepare_data(obs);
    CHECK(data.num_loci() == 4);
    CHECK(data.total_count() == 5.0);
    CHECK(data.num_pairs() == 4);
    // Sorted by (coverage, derived, frequency).
    CHECK(data.pair_coverage == std::vector<std::int64_t>({2, 3, 3, 4}));
    CHECK(data.pair_derived == std::vector<std::int64_t>({2, 1, 3, 4}));
    CHECK(data.counts == std::vector<double>({1, 2, 1, 1}));
    CHECK(data.frequencies == std::vector<double>({0.5, 0.25, 0.75, 0.875}));

    // Pairs are shared between loci with different frequencies.
    Observations obs2;
    obs2.push_back(5, 2, 0.1, 3);
    obs2.push_back(5, 2, 0.2, 4);
    obs2.push_back(5, 3, 0.2, 1);
    const PreparedData data2 = prepare_data(obs2);
    CHECK(data2.num_pairs() == 2);
    CHECK(data2.pair == std::vector<std::int32_t>({0, 0, 1}));

    Observations bad;
    bad.push_back(3, 4, 0.5, 1);
    CHECK_THROWS(prepare_data(bad), std::invalid_argument);
    Observations allfixed;
    allfixed.push_back(3, 1, 1.0, 1);
    allfixed.push_back(3, 1, 0.0, 1);
    CHECK_THROWS(prepare_data(allfixed), std::invalid_argument);
    Observations zerocoverage;
    zerocoverage.push_back(0, 0, 0.5, 1);
    CHECK_THROWS(prepare_data(zerocoverage), std::invalid_argument);
    Observations badfreq;
    badfreq.push_back(3, 1, 1.5, 1);
    CHECK_THROWS(prepare_data(badfreq), std::invalid_argument);
    CHECK_THROWS(prepare_data(Observations{}), std::invalid_argument);
}

void test_readers() {
    const fs::path dir = temp_dir();
    const std::string rows = "3,4,0.25,10\n0,2,0.5,1\n5,0,0.125,2\n";
    const std::vector<std::pair<std::string, std::string>> variants = {
        {"a.csv", "ancestral,derived,freq,counts\n" + rows},
        {"b.csv", rows},
        {"c.csv.gz", "ancestral,derived,freq,counts\n" + rows},
        {"d.txt", "3\t4\t0.25\t10\r\n0\t2\t0.5\t1\r\n\r\n5\t0\t0.125\t2\r\n"},
        {"e.txt.gz", "3 4   0.25 10\n0 2 0.5 1\n5 0 0.125 2"},  // no final newline
    };
    for (const auto& [name, text] : variants) {
        write_text(dir / name, text);
        ReadSummary summary;
        const Observations obs = read_dice_file((dir / name).string(), std::nullopt, &summary);
        CHECK(obs.coverages == std::vector<std::int64_t>({7, 2, 5}));
        CHECK(obs.derivedreads == std::vector<std::int64_t>({4, 2, 0}));
        CHECK(obs.frequencies == std::vector<double>({0.25, 0.5, 0.125}));
        CHECK(obs.counts == std::vector<std::int64_t>({10, 1, 2}));
        CHECK(summary.rows == 3 && summary.loci == 13);
    }
    // Explicit header settings.
    CHECK(read_dice_file((dir / "a.csv").string(), true).size() == 3);
    CHECK(read_dice_file((dir / "b.csv").string(), false).size() == 3);
    CHECK_THROWS(read_dice_file((dir / "a.csv").string(), false), std::runtime_error);  // header read as data

    // Base count and frequency files; identical rows are aggregated.
    write_text(dir / "bc.txt.gz", "4 7\n2 2\n4 7\n0 5\n");
    write_text(dir / "fr.txt", "freq\n0.25\n0.5\n0.25\n0.125\n");
    ReadSummary summary;
    const Observations obs =
        read_basecount_and_frequency_files((dir / "bc.txt.gz").string(), (dir / "fr.txt").string(), std::nullopt,
                                           std::nullopt, &summary);
    CHECK(summary.rows == 4 && !summary.header_first && summary.header_second);
    CHECK(obs.coverages == std::vector<std::int64_t>({7, 2, 5}));
    CHECK(obs.derivedreads == std::vector<std::int64_t>({4, 2, 0}));
    CHECK(obs.frequencies == std::vector<double>({0.25, 0.5, 0.125}));
    CHECK(obs.counts == std::vector<std::int64_t>({2, 1, 1}));

    // The bundled fixtures read the same way as through prepare_data.
    const PreparedData fixture = prepare_data(read_basecount_and_frequency_files(
        data_file("small_basecount.txt.gz"), data_file("small_freq.txt.gz"), std::nullopt, std::nullopt));
    CHECK(fixture.total_count() > 2000);

    // A UTF-8 byte order mark does not hide the first row.
    write_text(dir / "bom.csv", "\xEF\xBB\xBF" + rows);
    CHECK(read_dice_file((dir / "bom.csv").string(), std::nullopt).size() == 3);

    // As in Charon.jl, base count rows with zero coverage or an anchor frequency outside
    // (0, 1), including NaN, are dropped before the remaining rows are checked.
    write_text(dir / "bcjunk.txt", "4 7\n3 0\n1 5\n9 9\n");
    write_text(dir / "frjunk.txt", "0.25\n0.5\nNaN\n1.0\n");
    ReadSummary junk;
    const Observations kept = read_basecount_and_frequency_files(
        (dir / "bcjunk.txt").string(), (dir / "frjunk.txt").string(), std::nullopt, std::nullopt, &junk);
    CHECK(junk.rows == 4 && kept.size() == 1 && kept.coverages[0] == 7);
    write_text(dir / "frnone.txt", "0.0\n0.5\n1.0\n1.0\n");
    CHECK_THROWS(read_basecount_and_frequency_files((dir / "bcjunk.txt").string(), (dir / "frnone.txt").string(),
                                                    std::nullopt, std::nullopt),
                 std::runtime_error);

    // Errors.
    write_text(dir / "short.txt", "0.25\n0.5\n");
    CHECK_THROWS(read_basecount_and_frequency_files((dir / "bc.txt.gz").string(), (dir / "short.txt").string(),
                                                    std::nullopt, std::nullopt),
                 std::runtime_error);
    write_text(dir / "cols.csv", "3,4,0.25\n");
    CHECK_THROWS(read_dice_file((dir / "cols.csv").string(), std::nullopt), std::runtime_error);
    write_text(dir / "num.csv", "3,4,0.25,10\n3,x,0.25,10\n");
    CHECK_THROWS(read_dice_file((dir / "num.csv").string(), std::nullopt), std::runtime_error);
    write_text(dir / "range.csv", "3,4,1.25,10\n");
    CHECK_THROWS(read_dice_file((dir / "range.csv").string(), std::nullopt), std::runtime_error);
    write_text(dir / "derived.txt", "8 7\n");
    write_text(dir / "onefreq.txt", "0.5\n");
    CHECK_THROWS(read_basecount_and_frequency_files((dir / "derived.txt").string(), (dir / "onefreq.txt").string(),
                                                    std::nullopt, std::nullopt),
                 std::runtime_error);
    write_text(dir / "empty.csv", "\n\n");
    CHECK_THROWS(read_dice_file((dir / "empty.csv").string(), std::nullopt), std::runtime_error);
    CHECK_THROWS(read_dice_file((dir / "does_not_exist.csv").string(), std::nullopt), std::runtime_error);
    // A truncated gzip file is an error, not silently shorter data.
    {
        std::string big;
        for (int i = 0; i < 20000; ++i) big += "3,4,0.25," + std::to_string(i + 1) + "\n";
        write_text(dir / "full.csv.gz", big);
        std::ifstream in(dir / "full.csv.gz", std::ios::binary);
        std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        std::ofstream(dir / "truncated.csv.gz", std::ios::binary) << bytes.substr(0, bytes.size() / 2);
        CHECK_THROWS(read_dice_file((dir / "truncated.csv.gz").string(), std::nullopt), std::runtime_error);
    }
}

void test_format_number() {
    CHECK(format_number(0.25) == "0.25");
    CHECK(format_number(1.0) == "1.0");
    CHECK(format_number(-3.0) == "-3.0");
    CHECK(format_number(0.0) == "0.0");
    CHECK(format_number(1e-5) == "1.0e-5");
    CHECK(format_number(2.5e-7) == "2.5e-7");
    CHECK(format_number(1e20) == "1.0e20");
    CHECK(format_number(0.1 + 0.2) == "0.30000000000000004");
    CHECK(format_number(123456.0) == "123456.0");
    CHECK(format_number(kNegInf) == "-Inf");
    CHECK(format_number(std::numeric_limits<double>::quiet_NaN()) == "NaN");
    Rng rng(3);
    for (int i = 0; i < 10000; ++i) {
        const double x = (rng.uniform() - 0.5) * std::pow(10.0, static_cast<int>(rng.uniform() * 40) - 20);
        CHECK(std::strtod(format_number(x).c_str(), nullptr) == x);
    }
}

void test_priors() {
    const ContinuousPrior u = ContinuousPrior::parse("uniform:0:0.5");
    CHECK_CLOSE(u.logpdf(0.2), std::log(2.0), 1e-15);
    CHECK(u.logpdf(0.0) == std::log(2.0) && u.logpdf(0.5) == std::log(2.0));
    CHECK(u.logpdf(-0.01) == kNegInf && u.logpdf(0.51) == kNegInf);
    const ContinuousPrior b = ContinuousPrior::parse("beta:2:3");
    CHECK_CLOSE(b.logpdf(0.3), std::log(12 * 0.3 * 0.7 * 0.7), 1e-14);
    CHECK(b.logpdf(1.2) == kNegInf);
    const ContinuousPrior bs = ContinuousPrior::parse("beta:2:3:0:0.5");
    CHECK_CLOSE(bs.logpdf(0.15), std::log(12 * 0.3 * 0.7 * 0.7) + std::log(2.0), 1e-14);
    CHECK(bs.lower() == 0.0 && bs.upper() == 0.5);
    CHECK_CLOSE(ContinuousPrior::parse("beta:1:1").logpdf(0.0), 0.0, 1e-15);
    CHECK_THROWS(ContinuousPrior::parse("normal:0:1"), std::invalid_argument);
    CHECK_THROWS(ContinuousPrior::parse("uniform:1:0"), std::invalid_argument);
    CHECK_THROWS(ContinuousPrior::parse("uniform:0"), std::invalid_argument);
    CHECK_THROWS(ContinuousPrior::parse("beta:0:1"), std::invalid_argument);

    const DiscretePrior n = DiscretePrior::parse("uniform:1:10");
    CHECK(n.minimum() == 1 && n.maximum() == 10);
    CHECK_CLOSE(n.logpdf(4), -std::log(10.0), 1e-15);
    CHECK(n.logpdf(0) == kNegInf && n.logpdf(11) == kNegInf);
    const DiscretePrior n3 = DiscretePrior::parse("uniform:3:5");
    CHECK(n3.minimum() == 3 && n3.logpdf(2) == kNegInf && n3.logpdf(3) == -std::log(3.0));
    const DiscretePrior w = DiscretePrior::parse("weights:1,0,3");
    CHECK(w.minimum() == 1 && w.maximum() == 3);
    CHECK_CLOSE(w.logpdf(1), std::log(0.25), 1e-15);
    CHECK(w.logpdf(2) == kNegInf);
    CHECK_CLOSE(w.logpdf(3), std::log(0.75), 1e-15);
    CHECK_THROWS(DiscretePrior::parse("uniform:0:10"), std::invalid_argument);
    CHECK_THROWS(DiscretePrior::parse("weights:0,0"), std::invalid_argument);
    CHECK_THROWS(DiscretePrior::parse("weights:1,-1"), std::invalid_argument);
    CHECK_THROWS(DiscretePrior::parse("poisson:3"), std::invalid_argument);
}

void test_exactposterior() {
    const PreparedData data = small_data();
    Logger log(true);
    const std::vector<int> ns = {1, 3};
    const std::vector<double> tC = {0.1, 0.4}, tA = {0.0, 0.2, 0.5}, eps = {0.0, 0.05};
    const GridResult g = exactposterior(data, ns, tC, tA, eps, 2, 0, log);
    CHECK(g.loglik.size() == 24);
    std::size_t i = 0;
    for (int n : ns) {
        const ConditionalModel model(n, data);
        auto ws = model.make_workspace();
        for (double c : tC) {
            for (double a : tA) {
                for (double e : eps) {
                    CHECK(g.n[i] == n && g.tauC[i] == c && g.tauA[i] == a && g.eps[i] == e);
                    CHECK_CLOSE(g.loglik[i], model.loglik(c, a, e, ws), 1e-13);
                    ++i;
                }
            }
        }
    }
    // Thread count does not change the result.
    const GridResult g1 = exactposterior(data, ns, tC, tA, eps, 1, 0, log);
    CHECK(g1.loglik == g.loglik);
    CHECK_THROWS(exactposterior(data, {}, tC, tA, eps, 1, 0, log), std::invalid_argument);
    CHECK_THROWS(exactposterior(data, {0}, tC, tA, eps, 1, 0, log), std::invalid_argument);
    CHECK_THROWS(exactposterior(data, ns, {-0.1}, tA, eps, 1, 0, log), std::invalid_argument);
    CHECK_THROWS(exactposterior(data, ns, tC, tA, {0.5}, 1, 0, log), std::invalid_argument);
    CHECK_THROWS(exactposterior(data, {kMaxIndividualsGrid + 1}, tC, tA, eps, 1, 0, log), std::invalid_argument);
    const std::vector<double> huge(100000, 0.5);
    CHECK_THROWS(exactposterior(data, ns, huge, huge, huge, 1, 0, log), std::invalid_argument);
    // The largest allowed n still gives finite, accurate values.
    const GridResult big = exactposterior(data, {kMaxIndividualsGrid}, {0.05}, {0.1}, {0.01}, 1, 0, log);
    CHECK(std::isfinite(big.loglik[0]));
}

void check_chain_invariants(const Chain& ch, const MCMCSettings& s) {
    for (int step = 0; step < ch.nsteps; ++step) {
        const int n = ch.nsample[step];
        CHECK(n >= 1 && n <= ch.nmax);
        CHECK(s.prioronn.logpdf(n) > kNegInf);
        for (int k = 1; k <= ch.nmax; ++k) {
            const std::size_t idx = static_cast<std::size_t>(step) * ch.nmax + (k - 1);
            CHECK(s.prioronTauC.logpdf(ch.tauC[idx]) > kNegInf);
            CHECK(s.prioronTauA.logpdf(ch.tauA[idx]) > kNegInf);
            CHECK(s.prioronEps.logpdf(ch.eps[idx]) > kNegInf);
            if (step > 0 && !ch.accepted[static_cast<std::size_t>(step) * (ch.nmax + 1) + k]) {
                CHECK(ch.tauC[idx] == ch.tauC[idx - ch.nmax]);
            }
        }
        const double lj = ch.logjointprob[static_cast<std::size_t>(step) * (ch.nmax + 1)];
        CHECK(std::isfinite(lj));
        CHECK_CLOSE(lj,
                    ch.logjointprob[static_cast<std::size_t>(step) * (ch.nmax + 1) + n] + s.prioronn.logpdf(n),
                    1e-12);
    }
}

void test_mcmc() {
    const PreparedData data = small_data();
    Logger log(true);
    MCMCSettings s;
    s.nchains = 3;
    s.nsteps = 400;
    s.prioronn = DiscretePrior::uniform(1, 4);
    s.seed = 42;
    s.messages = 0;
    s.threads = 1;
    const std::vector<Chain> a = mcmc_sampler(data, s, log);
    s.threads = 3;
    const std::vector<Chain> b = mcmc_sampler(data, s, log);
    CHECK(a.size() == 3 && b.size() == 3);
    for (int c = 0; c < 3; ++c) {
        // Reproducible for a seed, whatever the number of threads.
        CHECK(a[c].nsample == b[c].nsample && a[c].tauC == b[c].tauC && a[c].eps == b[c].eps &&
              a[c].logjointprob == b[c].logjointprob);
        check_chain_invariants(a[c], s);
    }
    CHECK(a[0].tauC != a[1].tauC);  // different chains get different streams
    // Spare threads evaluate the conditional likelihoods of a chain in parallel; the
    // chain must not depend on that.
    s.nchains = 1;
    s.threads = 1;
    const Chain serial = mcmc_sampler(data, s, log)[0];
    s.threads = 4;
    const Chain team = mcmc_sampler(data, s, log)[0];
    CHECK(serial.nsample == team.nsample && serial.tauC == team.tauC && serial.tauA == team.tauA &&
          serial.eps == team.eps && serial.accepted == team.accepted && serial.logjointprob == team.logjointprob);
    CHECK(serial.tauC == a[0].tauC);
    s.nchains = 3;
    s.seed = 43;
    CHECK(mcmc_sampler(data, s, log)[0].tauC != a[0].tauC);

    // Priors on n with a lower bound > 1 and with zero weights (Charon.jl fails on these).
    s.nchains = 1;
    s.nsteps = 300;
    s.prioronn = DiscretePrior::uniform(3, 5);
    check_chain_invariants(mcmc_sampler(data, s, log)[0], s);
    s.prioronn = DiscretePrior::weights({1, 0, 1});
    check_chain_invariants(mcmc_sampler(data, s, log)[0], s);
    s.prioronn = DiscretePrior::uniform(1, 1);  // a single n
    const Chain single = mcmc_sampler(data, s, log)[0];
    check_chain_invariants(single, s);

    // A beta prior with B < 1 has infinite density at its upper bound, which is on the grid
    // for the starting values; chains must not start (and get stuck) there.
    s.prioronn = DiscretePrior::uniform(1, 3);
    s.nsteps = 400;
    s.prioronTauA = ContinuousPrior::beta(0.5, 0.5, 0.0, 0.3);
    s.prioronTauC = ContinuousPrior::beta(2.0, 0.5, 0.0, 0.5);
    const Chain betachain = mcmc_sampler(data, s, log)[0];
    check_chain_invariants(betachain, s);
    int moved = 0;
    for (int step = 1; step < s.nsteps; ++step) moved += betachain.accepted[static_cast<std::size_t>(step) * 4 + 1];
    CHECK(moved > 0);
    s.prioronTauA = ContinuousPrior::uniform(0.0, 1.0);
    s.prioronTauC = ContinuousPrior::uniform(0.0, 1.0);

    // Invalid settings.
    s.prioronn = DiscretePrior::uniform(1, 11);
    CHECK_THROWS(mcmc_sampler(data, s, log), std::invalid_argument);
    s.prioronn = DiscretePrior::uniform(1, 3);
    s.prioronEps = ContinuousPrior::uniform(0.0, 0.6);
    CHECK_THROWS(mcmc_sampler(data, s, log), std::invalid_argument);
    s.prioronEps = ContinuousPrior::uniform(0.0, 0.5);
    s.nsteps = 0;
    CHECK_THROWS(mcmc_sampler(data, s, log), std::invalid_argument);
}

void test_mcmc_adaptation() {
    // A long chain with a single n exercises the adaptation of the proposals,
    // including the scale factor updates (10000 < step < nsteps/4).
    const PreparedData data = small_data();
    Logger log(true);
    MCMCSettings s;
    s.nchains = 1;
    s.nsteps = 48000;
    s.prioronn = DiscretePrior::uniform(1, 1);
    s.seed = 5;
    s.messages = 0;
    s.scalingmessages = true;
    const Chain ch = mcmc_sampler(data, s, log)[0];
    check_chain_invariants(ch, s);
    int accepted = 0;
    for (int step = 24000; step < 48000; ++step) accepted += ch.accepted[static_cast<std::size_t>(step) * 2 + 1];
    const double rate = accepted / 24000.0;
    CHECK(rate > 0.03 && rate < 0.6);
}

void test_unpackposterior() {
    // Port of the unpackposterior test.
    const int nsteps = 100, nmax = 10;
    Chain ch(nsteps, nmax);
    for (int step = 0; step < nsteps; ++step) {
        ch.nsample[step] = step % 10 + 1;  // repeat(1:10, 10)
        for (int n = 1; n <= nmax; ++n) {
            const double v = step * nmax + n;
            ch.tauC_at(step, n) = v;
            ch.tauA_at(step, n) = -v;
            ch.eps_at(step, n) = v * v;
        }
        for (int col = 0; col <= nmax; ++col) ch.logjointprob_at(step, col) = -step - 0.5 * col;
    }
    const fs::path path = temp_dir() / "unpacked.csv";
    {
        OutputFile out(path.string());
        write_unpacked_posterior({ch, ch}, out);
        out.close();
    }
    const auto rows = read_csv(path);
    CHECK(rows.size() == 1 + 2 * nsteps);
    CHECK(rows[0] == std::vector<std::string>({"nsample", "τCsample", "τAsample", "ϵsample", "accepted",
                                                "logjointprob", "chainid"}));
    for (int i = 0; i < 2 * nsteps; ++i) {
        const auto& r = rows[1 + i];
        const int step = i % nsteps;
        const int n = std::stoi(r[0]);
        CHECK(n == step % 10 + 1);
        CHECK(std::stod(r[1]) == step * nmax + n);
        CHECK(std::stod(r[2]) == -(step * nmax + n));
        CHECK(std::stod(r[3]) == std::stod(r[1]) * std::stod(r[1]));
        CHECK(r[4] == "true");
        CHECK(std::stod(r[5]) == -step);
        CHECK(r[6] == (i < nsteps ? "1" : "2"));
    }
    // accepted is false exactly when nothing changed.
    Chain still(3, 2);
    for (int step = 0; step < 3; ++step) {
        still.nsample[step] = 2;
        for (int n = 1; n <= 2; ++n) still.tauC_at(step, n) = still.tauA_at(step, n) = still.eps_at(step, n) = 0.1;
    }
    still.tauA_at(2, 2) = 0.2;
    const fs::path path2 = temp_dir() / "still.csv.gz";
    {
        OutputFile out(path2.string());
        write_unpacked_posterior({still}, out);
        out.close();
    }
    const auto rows2 = read_csv(path2);  // also checks that .gz output is readable
    CHECK(rows2.size() == 4);
    CHECK(rows2[1][4] == "true" && rows2[2][4] == "false" && rows2[3][4] == "true");
}

void test_rng() {
    Rng a(1), b(1), c(2);
    bool differs = false;
    for (int i = 0; i < 100; ++i) {
        const std::uint64_t x = a.next();
        CHECK(x == b.next());
        differs = differs || x != c.next();
    }
    CHECK(differs);
    Rng r(9);
    double sum = 0, sumsq = 0, umin = 1, umax = 0;
    const int N = 200000;
    for (int i = 0; i < N; ++i) {
        const double z = r.normal();
        sum += z;
        sumsq += z * z;
        const double u = r.uniform();
        umin = std::min(umin, u);
        umax = std::max(umax, u);
    }
    CHECK(std::fabs(sum / N) < 0.01);
    CHECK(std::fabs(sumsq / N - 1.0) < 0.02);
    CHECK(umin >= 0.0 && umax < 1.0);
}

}  // namespace

int main() {
    const std::vector<std::pair<const char*, std::function<void()>>> tests = {
        {"makeQ", test_makeQ},
        {"EigenExpansion", test_eigen_expansion},
        {"binomial", test_binomial},
        {"loglik reference values (Charon.jl)", test_loglik_reference_values},
        {"loglik without drift", test_loglik_without_drift},
        {"prepare_data", test_prepare_data},
        {"readers", test_readers},
        {"format_number", test_format_number},
        {"priors", test_priors},
        {"exactposterior", test_exactposterior},
        {"mcmc", test_mcmc},
        {"mcmc adaptation", test_mcmc_adaptation},
        {"unpackposterior", test_unpackposterior},
        {"rng", test_rng},
    };
    for (const auto& [name, test] : tests) {
        const int before = g_failures;
        try {
            test();
        } catch (const std::exception& e) {
            ++g_failures;
            std::cerr << "unexpected exception in " << name << ": " << e.what() << "\n";
        }
        std::cout << (g_failures == before ? "[ OK ] " : "[FAIL] ") << name << "\n";
    }
    std::error_code ec;
    fs::remove_all(temp_dir(), ec);
    std::cout << g_checks << " checks, " << g_failures << " failures\n";
    return g_failures == 0 ? 0 : 1;
}
