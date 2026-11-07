#pragma once

// Priors on n, tauC, tauA and eps. Charon.jl accepts any Distributions.jl
// distribution; the command line version supports the distributions below, which
// are specified as strings such as "uniform:0:1" or "beta:2:5:0:0.5".

#include <string>
#include <vector>

namespace charon {

/// A continuous prior on a bounded interval [lower, upper]:
///  * "uniform:LO:HI"         Uniform(LO, HI)
///  * "beta:A:B"              Beta(A, B) on [0, 1]
///  * "beta:A:B:LO:HI"        Beta(A, B) rescaled to [LO, HI]
class ContinuousPrior {
public:
    static ContinuousPrior uniform(double lower, double upper);
    static ContinuousPrior beta(double a, double b, double lower = 0.0, double upper = 1.0);
    /// Throws std::invalid_argument on a malformed specification.
    static ContinuousPrior parse(const std::string& spec);

    double logpdf(double x) const;
    double lower() const { return lower_; }
    double upper() const { return upper_; }
    std::string describe() const;

private:
    enum class Kind { Uniform, Beta };
    ContinuousPrior(Kind kind, double lower, double upper, double a, double b);

    Kind kind_;
    double lower_;
    double upper_;
    double a_;
    double b_;
    double lognorm_;  // log normalising constant
};

/// A prior on the number of individuals n on {1, 2, ..., maximum()}:
///  * "uniform:A:B"   discrete uniform on {A, ..., B}
///  * "weights:w1,w2,...,wK"   P(n = k) proportional to wk (like Categorical in Julia)
class DiscretePrior {
public:
    static DiscretePrior uniform(int lower, int upper);
    static DiscretePrior weights(const std::vector<double>& weights);
    static DiscretePrior parse(const std::string& spec);

    double logpdf(int n) const;
    /// Smallest value of the support as reported by Distributions.jl's minimum().
    int minimum() const { return minimum_; }
    int maximum() const { return maximum_; }
    std::string describe() const;

private:
    DiscretePrior() = default;

    std::vector<double> logp_;  // logp_[n] for n = 0, ..., maximum_
    int minimum_ = 1;
    int maximum_ = 1;
    std::string description_;
};

/// Split text at separator and parse every part as a number.
std::vector<double> split_numbers(const std::string& text, char separator);

}  // namespace charon
