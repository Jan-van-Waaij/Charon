#pragma once

// The grid posterior (exactposterior), the MCMC sampler (MCMCsampler) and the
// unconditional posterior (unpackposterior) of Charon.jl.

#include <cstdint>
#include <vector>

#include "charon/io.hpp"
#include "charon/model.hpp"
#include "charon/priors.hpp"

namespace charon {

/// Largest number of individuals supported by the MCMC sampler (as in Charon.jl).
inline constexpr int kMaxIndividualsMCMC = 10;
/// Largest number of individuals supported by exactposterior. The eigenvectors of Q and
/// Qꜜ become too ill-conditioned for accurate matrix exponentials beyond this.
inline constexpr int kMaxIndividualsGrid = 10;
/// Largest number of grid points for exactposterior.
inline constexpr std::size_t kMaxGridPoints = 200000000;

/// Number of worker threads to use when the user asked for `requested` (0 = all cores).
int resolve_threads(int requested);

/// The log-likelihood for every combination (n, tauC, tauA, eps) of the grid, in the
/// order of Charon.jl: n outermost, then tauC, then tauA, eps innermost.
struct GridResult {
    std::vector<int> n;
    std::vector<double> tauC;
    std::vector<double> tauA;
    std::vector<double> eps;
    std::vector<double> loglik;
};

/// Port of exactposterior. messages > 0 prints progress roughly every `messages`
/// evaluations. Throws std::invalid_argument on invalid grids.
GridResult exactposterior(const PreparedData& data, const std::vector<int>& nrange, const std::vector<double>& tauCrange,
                          const std::vector<double>& tauArange, const std::vector<double>& epsrange, int threads,
                          std::int64_t messages, Logger& log);

struct MCMCSettings {
    int nchains = 4;
    int nsteps = 100000;
    DiscretePrior prioronn = DiscretePrior::uniform(1, 10);
    ContinuousPrior prioronTauC = ContinuousPrior::uniform(0.0, 1.0);
    ContinuousPrior prioronTauA = ContinuousPrior::uniform(0.0, 1.0);
    ContinuousPrior prioronEps = ContinuousPrior::uniform(0.0, 0.5);
    std::int64_t messages = -1;   ///< progress every `messages` steps; negative: nsteps / 100; 0: none
    bool scalingmessages = false; ///< report changes of the adaptive scale factors
    int threads = 0;              ///< 0: all cores
    std::uint64_t seed = 0;
};

/// One Markov chain, as returned by Charon.jl's MCMCsampler. Matrices are stored
/// row-major with one row per step; column n-1 belongs to the chain conditioned on n.
struct Chain {
    int nsteps = 0;
    int nmax = 0;
    std::vector<int> nsample;            ///< nsteps
    std::vector<double> tauC;            ///< nsteps x nmax
    std::vector<double> tauA;            ///< nsteps x nmax
    std::vector<double> eps;             ///< nsteps x nmax
    std::vector<std::uint8_t> accepted;  ///< nsteps x (nmax+1); column 0 is the n chain
    std::vector<double> logjointprob;    ///< nsteps x (nmax+1); column 0 is the unconditional chain

    Chain(int nsteps, int nmax);
    double& tauC_at(int step, int n) { return tauC[static_cast<std::size_t>(step) * nmax + (n - 1)]; }
    double& tauA_at(int step, int n) { return tauA[static_cast<std::size_t>(step) * nmax + (n - 1)]; }
    double& eps_at(int step, int n) { return eps[static_cast<std::size_t>(step) * nmax + (n - 1)]; }
    std::uint8_t& accepted_at(int step, int col) { return accepted[static_cast<std::size_t>(step) * (nmax + 1) + col]; }
    double& logjointprob_at(int step, int col) {
        return logjointprob[static_cast<std::size_t>(step) * (nmax + 1) + col];
    }
};

/// Port of MCMCsampler: validates the settings, finds starting values on a grid and
/// runs settings.nchains chains in parallel. Throws std::invalid_argument on invalid
/// settings and std::runtime_error when the data do not support the initial values.
std::vector<Chain> mcmc_sampler(const PreparedData& data, const MCMCSettings& settings, Logger& log);

/// Port of unpackposterior: writes the unconditional posterior as CSV with columns
/// nsample, τCsample, τAsample, ϵsample, accepted, logjointprob, chainid.
void write_unpacked_posterior(const std::vector<Chain>& chains, OutputFile& out);

/// Writes the full chain, including the chains conditioned on every n.
void write_raw_chain(const Chain& chain, OutputFile& out);

/// Writes the result of exactposterior as CSV with columns n, τC, τA, ϵ, loglik.
void write_grid(const GridResult& grid, OutputFile& out);

}  // namespace charon
