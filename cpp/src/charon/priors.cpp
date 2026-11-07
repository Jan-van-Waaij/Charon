#include "charon/priors.hpp"

#include <cmath>
#include <cstdlib>
#include <limits>
#include <sstream>
#include <stdexcept>

#include "charon/io.hpp"

namespace charon {

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

double parse_number(const std::string& token, const std::string& spec) {
    char* end = nullptr;
    const double v = std::strtod(token.c_str(), &end);
    if (token.empty() || end != token.c_str() + token.size() || !std::isfinite(v)) {
        throw std::invalid_argument("Could not read '" + token + "' as a number in '" + spec + "'.");
    }
    return v;
}

std::vector<std::string> split(const std::string& text, char separator) {
    std::vector<std::string> parts;
    std::string current;
    std::istringstream stream(text);
    while (std::getline(stream, current, separator)) parts.push_back(current);
    if (!text.empty() && text.back() == separator) parts.emplace_back();
    return parts;
}

}  // namespace

std::vector<double> split_numbers(const std::string& text, char separator) {
    std::vector<double> numbers;
    for (const std::string& part : split(text, separator)) numbers.push_back(parse_number(part, text));
    return numbers;
}

ContinuousPrior::ContinuousPrior(Kind kind, double lower, double upper, double a, double b)
    : kind_(kind), lower_(lower), upper_(upper), a_(a), b_(b), lognorm_(0.0) {
    if (!(std::isfinite(lower) && std::isfinite(upper) && lower < upper)) {
        throw std::invalid_argument("A prior needs finite bounds with lower < upper.");
    }
    if (kind == Kind::Uniform) {
        lognorm_ = -std::log(upper - lower);
    } else {
        if (!(a > 0.0 && b > 0.0 && std::isfinite(a) && std::isfinite(b))) {
            throw std::invalid_argument("The parameters of a beta prior should be positive.");
        }
        lognorm_ = std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) - std::log(upper - lower);
    }
}

ContinuousPrior ContinuousPrior::uniform(double lower, double upper) {
    return ContinuousPrior(Kind::Uniform, lower, upper, 1.0, 1.0);
}

ContinuousPrior ContinuousPrior::beta(double a, double b, double lower, double upper) {
    return ContinuousPrior(Kind::Beta, lower, upper, a, b);
}

ContinuousPrior ContinuousPrior::parse(const std::string& spec) {
    const std::vector<std::string> parts = split(spec, ':');
    if (parts.empty()) throw std::invalid_argument("Empty prior specification.");
    const std::string& name = parts[0];
    if (name == "uniform" && parts.size() == 3) {
        return uniform(parse_number(parts[1], spec), parse_number(parts[2], spec));
    }
    if (name == "beta" && parts.size() == 3) {
        return beta(parse_number(parts[1], spec), parse_number(parts[2], spec));
    }
    if (name == "beta" && parts.size() == 5) {
        return beta(parse_number(parts[1], spec), parse_number(parts[2], spec), parse_number(parts[3], spec),
                    parse_number(parts[4], spec));
    }
    throw std::invalid_argument("Unknown prior '" + spec +
                                "'. Use uniform:LO:HI, beta:A:B or beta:A:B:LO:HI.");
}

double ContinuousPrior::logpdf(double x) const {
    if (!(x >= lower_ && x <= upper_)) return kNegInf;
    if (kind_ == Kind::Uniform) return lognorm_;
    const double z = (x - lower_) / (upper_ - lower_);
    // Treat (a-1) log z as 0 when a == 1, also at z == 0 (as Distributions.jl does).
    const double left = a_ == 1.0 ? 0.0 : (a_ - 1.0) * std::log(z);
    const double right = b_ == 1.0 ? 0.0 : (b_ - 1.0) * std::log1p(-z);
    return lognorm_ + left + right;
}

std::string ContinuousPrior::describe() const {
    if (kind_ == Kind::Uniform) return "Uniform(" + format_number(lower_) + ", " + format_number(upper_) + ")";
    std::string s = "Beta(" + format_number(a_) + ", " + format_number(b_) + ")";
    if (lower_ != 0.0 || upper_ != 1.0) s += " on [" + format_number(lower_) + ", " + format_number(upper_) + "]";
    return s;
}

DiscretePrior DiscretePrior::uniform(int lower, int upper) {
    if (lower > upper) throw std::invalid_argument("The prior on n needs lower <= upper.");
    if (lower < 1) {
        throw std::invalid_argument(
            "The prior on n should have support on positive integers. There is always at least one individual.");
    }
    DiscretePrior prior;
    prior.minimum_ = lower;
    prior.maximum_ = upper;
    prior.logp_.assign(static_cast<std::size_t>(upper) + 1, kNegInf);
    const double lp = -std::log(static_cast<double>(upper - lower + 1));
    for (int n = lower; n <= upper; ++n) prior.logp_[n] = lp;
    prior.description_ = "DiscreteUniform(" + std::to_string(lower) + ", " + std::to_string(upper) + ")";
    return prior;
}

DiscretePrior DiscretePrior::weights(const std::vector<double>& weights) {
    if (weights.empty()) throw std::invalid_argument("The prior on n needs at least one weight.");
    double total = 0.0;
    for (double w : weights) {
        if (!(w >= 0.0) || !std::isfinite(w)) throw std::invalid_argument("The weights of the prior on n should be non-negative.");
        total += w;
    }
    if (!(total > 0.0)) throw std::invalid_argument("At least one weight of the prior on n should be positive.");
    DiscretePrior prior;
    prior.minimum_ = 1;  // like Categorical in Distributions.jl
    prior.maximum_ = static_cast<int>(weights.size());
    prior.logp_.assign(weights.size() + 1, kNegInf);
    std::string description = "Categorical([";
    for (std::size_t k = 0; k < weights.size(); ++k) {
        prior.logp_[k + 1] = weights[k] > 0.0 ? std::log(weights[k] / total) : kNegInf;
        if (k > 0) description += ", ";
        description += format_number(weights[k] / total);
    }
    prior.description_ = description + "])";
    return prior;
}

DiscretePrior DiscretePrior::parse(const std::string& spec) {
    const std::size_t colon = spec.find(':');
    const std::string name = spec.substr(0, colon);
    const std::string rest = colon == std::string::npos ? "" : spec.substr(colon + 1);
    if (name == "uniform") {
        const std::vector<double> bounds = split_numbers(rest, ':');
        if (bounds.size() == 2 && bounds[0] == std::floor(bounds[0]) && bounds[1] == std::floor(bounds[1]) &&
            std::fabs(bounds[0]) < 1e6 && std::fabs(bounds[1]) < 1e6) {
            return uniform(static_cast<int>(bounds[0]), static_cast<int>(bounds[1]));
        }
    } else if (name == "weights") {
        return weights(split_numbers(rest, ','));
    }
    throw std::invalid_argument("Unknown prior on n '" + spec + "'. Use uniform:A:B or weights:w1,w2,...");
}

double DiscretePrior::logpdf(int n) const {
    if (n < 0 || n > maximum_) return kNegInf;
    return logp_[static_cast<std::size_t>(n)];
}

std::string DiscretePrior::describe() const { return description_; }

}  // namespace charon
