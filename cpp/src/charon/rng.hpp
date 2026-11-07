#pragma once

#include <cmath>
#include <cstdint>

namespace charon {

/// splitmix64, used to expand a 64-bit seed into generator state.
inline std::uint64_t splitmix64(std::uint64_t& x) {
    std::uint64_t z = (x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

/// xoshiro256++ (the generator Julia uses by default), with its own uniform and
/// normal samplers so that a seed gives identical chains on every platform.
class Rng {
public:
    explicit Rng(std::uint64_t seed) {
        for (auto& word : s_) word = splitmix64(seed);
    }

    std::uint64_t next() {
        const std::uint64_t result = rotl(s_[0] + s_[3], 23) + s_[0];
        const std::uint64_t t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return result;
    }

    /// Uniform on [0, 1) with 53 random bits.
    double uniform() { return static_cast<double>(next() >> 11) * 0x1.0p-53; }

    /// Standard normal (Marsaglia polar method).
    double normal() {
        if (has_spare_) {
            has_spare_ = false;
            return spare_;
        }
        double u, v, s;
        do {
            u = 2.0 * uniform() - 1.0;
            v = 2.0 * uniform() - 1.0;
            s = u * u + v * v;
        } while (s >= 1.0 || s == 0.0);
        const double factor = std::sqrt(-2.0 * std::log(s) / s);
        spare_ = v * factor;
        has_spare_ = true;
        return u * factor;
    }

private:
    static std::uint64_t rotl(std::uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }

    std::uint64_t s_[4];
    double spare_ = 0.0;
    bool has_spare_ = false;
};

}  // namespace charon
