#include "charon/sampler.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <condition_variable>
#include <exception>
#include <functional>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>

#include "charon/rng.hpp"

namespace charon {

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

/// Runs body(0), ..., body(count-1) on up to `threads` threads; rethrows the first exception.
template <class Body>
void parallel_for(std::size_t count, int threads, Body&& body) {
    const int workers = static_cast<int>(std::max<std::size_t>(1, std::min<std::size_t>(threads, count)));
    if (workers == 1) {
        for (std::size_t i = 0; i < count; ++i) body(i);
        return;
    }
    std::atomic<std::size_t> next{0};
    std::atomic<bool> failed{false};
    std::exception_ptr error;
    std::mutex error_mutex;
    auto worker = [&]() {
        while (!failed.load()) {
            const std::size_t i = next.fetch_add(1);
            if (i >= count) return;
            try {
                body(i);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
                failed.store(true);
                return;
            }
        }
    };
    std::vector<std::thread> pool;
    pool.reserve(workers - 1);
    for (int t = 0; t < workers - 1; ++t) pool.emplace_back(worker);
    worker();
    for (std::thread& t : pool) t.join();
    if (error) std::rethrow_exception(error);
}

/// A fixed team of threads that repeatedly runs small batches of tasks, used to
/// evaluate the conditional likelihoods of one MCMC step in parallel. The calling
/// thread takes part in the work.
class WorkerTeam {
public:
    explicit WorkerTeam(int threads) {
        for (int t = 1; t < threads; ++t) helpers_.emplace_back([this] { helper_loop(); });
    }
    ~WorkerTeam() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
        }
        wake_.notify_all();
        for (std::thread& t : helpers_) t.join();
    }
    WorkerTeam(const WorkerTeam&) = delete;
    WorkerTeam& operator=(const WorkerTeam&) = delete;

    /// Runs body(0), ..., body(count-1) and returns when all are done. body must not throw.
    void run(std::size_t count, const std::function<void(std::size_t)>& body) {
        if (helpers_.empty() || count <= 1) {
            for (std::size_t i = 0; i < count; ++i) body(i);
            return;
        }
        {
            std::lock_guard<std::mutex> lock(mutex_);
            body_ = &body;
            count_ = count;
            next_.store(0);
            busy_.store(static_cast<int>(helpers_.size()));
            ++generation_;
        }
        wake_.notify_all();
        work();
        while (busy_.load(std::memory_order_acquire) > 0) std::this_thread::yield();
    }

private:
    void work() {
        for (std::size_t i = next_.fetch_add(1); i < count_; i = next_.fetch_add(1)) (*body_)(i);
    }

    void helper_loop() {
        std::uint64_t seen = 0;
        for (;;) {
            {
                std::unique_lock<std::mutex> lock(mutex_);
                wake_.wait(lock, [&] { return stop_ || generation_ != seen; });
                if (stop_) return;
                seen = generation_;
            }
            work();
            busy_.fetch_sub(1, std::memory_order_release);
        }
    }

    std::vector<std::thread> helpers_;
    std::mutex mutex_;
    std::condition_variable wake_;
    std::uint64_t generation_ = 0;
    bool stop_ = false;
    const std::function<void(std::size_t)>* body_ = nullptr;
    std::size_t count_ = 0;
    std::atomic<std::size_t> next_{0};
    std::atomic<int> busy_{0};
};

using ModelTable = std::vector<std::unique_ptr<ConditionalModel>>;  // indexed by n

ModelTable build_models(const PreparedData& data, const std::vector<int>& ns, int threads) {
    const int nmax = *std::max_element(ns.begin(), ns.end());
    ModelTable models(static_cast<std::size_t>(nmax) + 1);
    std::vector<int> distinct(ns);
    std::sort(distinct.begin(), distinct.end());
    distinct.erase(std::unique(distinct.begin(), distinct.end()), distinct.end());
    parallel_for(distinct.size(), threads, [&](std::size_t i) {
        const int n = distinct[i];
        models[n] = std::make_unique<ConditionalModel>(n, data);
    });
    return models;
}

/// The grid evaluation shared by exactposterior and the initialisation of the sampler.
GridResult evaluate_grid(const ModelTable& models, const std::vector<int>& nrange, const std::vector<double>& tauCrange,
                         const std::vector<double>& tauArange, const std::vector<double>& epsrange, int threads,
                         std::int64_t messages, Logger& log) {
    const std::size_t NN = nrange.size(), NC = tauCrange.size(), NA = tauArange.size(), NE = epsrange.size();
    const std::size_t nrows = NN * NC * NA * NE;
    GridResult g;
    g.n.resize(nrows);
    g.tauC.resize(nrows);
    g.tauA.resize(nrows);
    g.eps.resize(nrows);
    g.loglik.resize(nrows);
    std::atomic<std::size_t> done{0};
    parallel_for(NN * NC, threads, [&](std::size_t task) {
        const std::size_t a = task / NC;
        const std::size_t iC = task % NC;
        const ConditionalModel& model = *models[nrange[a]];
        ConditionalModel::Workspace ws = model.make_workspace();
        Mat expC, expA, M, tmp;
        model.exp_Qdown(tauCrange[iC], expC, tmp);
        for (std::size_t iA = 0; iA < NA; ++iA) {
            model.exp_Q(tauArange[iA], expA, tmp);
            M.noalias() = expA * expC;
            for (std::size_t iE = 0; iE < NE; ++iE) {
                const std::size_t index = ((a * NC + iC) * NA + iA) * NE + iE;
                g.n[index] = nrange[a];
                g.tauC[index] = tauCrange[iC];
                g.tauA[index] = tauArange[iA];
                g.eps[index] = epsrange[iE];
                g.loglik[index] = model.loglik_given_matrix(M, epsrange[iE], ws);
            }
        }
        if (messages > 0) {
            const std::size_t before = done.fetch_add(NA * NE);
            const std::size_t after = before + NA * NE;
            const auto every = static_cast<std::size_t>(messages);
            if (after / every > before / every) {
                log.info("I did " + std::to_string(after) + " posterior calculations, " +
                         std::to_string(100 * after / nrows) + "% finished!");
            }
        }
    });
    return g;
}

void check_grid(const std::vector<int>& nrange, const std::vector<double>& tauCrange,
                const std::vector<double>& tauArange, const std::vector<double>& epsrange) {
    if (nrange.empty()) throw std::invalid_argument("nrange should have at least one element.");
    if (*std::min_element(nrange.begin(), nrange.end()) <= 0) {
        throw std::invalid_argument("all n values in nrange should be positive.");
    }
    if (*std::max_element(nrange.begin(), nrange.end()) > kMaxIndividualsGrid) {
        throw std::invalid_argument("all n values in nrange should be at most " +
                                    std::to_string(kMaxIndividualsGrid) + ".");
    }
    auto check_tau = [](const std::vector<double>& range, const char* name) {
        if (range.empty()) throw std::invalid_argument(std::string(name) + " should have at least one element.");
        for (double t : range) {
            if (!(t >= 0.0) || !std::isfinite(t)) {
                throw std::invalid_argument(std::string("All values in ") + name + " should be non-negative.");
            }
        }
    };
    // Each grid point needs about 40 bytes of output; refuse grids that cannot fit.
    const long double points = static_cast<long double>(nrange.size()) * tauCrange.size() * tauArange.size() *
                               std::max<std::size_t>(epsrange.size(), 1);
    if (points > static_cast<long double>(kMaxGridPoints)) {
        throw std::invalid_argument("The grid has too many points (at most " + std::to_string(kMaxGridPoints) +
                                    " are allowed).");
    }
    check_tau(tauCrange, "τCrange");
    check_tau(tauArange, "τArange");
    if (epsrange.empty()) throw std::invalid_argument("ϵrange should have at least one element.");
    for (double e : epsrange) {
        if (!(e >= 0.0 && e < 0.5)) throw std::invalid_argument("All values in ϵrange should be in [0.0,0.5).");
    }
}

std::string format_vector(const std::vector<double>& v) {
    std::string s = "[";
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (i > 0) s += ", ";
        s += format_number(v[i]);
    }
    return s + "]";
}

/// Mean and corrected variance of column n of a chain matrix over rows first, first+stride, ..., <= last.
void strided_moments(const std::vector<double>& matrix, int ncols, int col, int first, int last, int stride,
                     double& mean, double& var) {
    double sum = 0.0;
    int count = 0;
    for (int r = first; r <= last; r += stride) {
        sum += matrix[static_cast<std::size_t>(r) * ncols + col];
        ++count;
    }
    mean = sum / count;
    double ss = 0.0;
    for (int r = first; r <= last; r += stride) {
        const double d = matrix[static_cast<std::size_t>(r) * ncols + col] - mean;
        ss += d * d;
    }
    var = count > 1 ? ss / (count - 1) : std::numeric_limits<double>::quiet_NaN();
}

struct StartingValues {
    int n = 1;
    std::vector<double> tauC, tauA, eps;  // index n-1
};

void run_chain(int chainindex, Chain& chain, const ModelTable& models, const StartingValues& start,
               const MCMCSettings& s, std::int64_t messages, std::uint64_t seed, int threads, Logger& log) {
    const int nmax = chain.nmax;
    const int nsteps = chain.nsteps;
    Rng rng(seed);
    // The likelihoods for n = 1, ..., nmax are independent; with spare threads they are
    // evaluated in parallel. All random numbers are drawn on this thread, so the chain
    // does not depend on the number of threads. Largest n (the most work) first.
    WorkerTeam team(std::min(threads, nmax));
    std::vector<ConditionalModel::Workspace> ws;
    ws.reserve(nmax);
    for (int n = 1; n <= nmax; ++n) ws.push_back(models[n]->make_workspace());
    auto logprior = [&](double tauC, double tauA, double eps) {
        return s.prioronTauC.logpdf(tauC) + s.prioronTauA.logpdf(tauA) + s.prioronEps.logpdf(eps);
    };

    // Initialise.
    int nprevious = start.n;
    std::vector<double> tauCprevious = start.tauC, tauAprevious = start.tauA, epsprevious = start.eps;
    chain.nsample[0] = nprevious;
    for (int n = 1; n <= nmax; ++n) {
        chain.tauC_at(0, n) = tauCprevious[n - 1];
        chain.tauA_at(0, n) = tauAprevious[n - 1];
        chain.eps_at(0, n) = epsprevious[n - 1];
    }
    for (int col = 0; col <= nmax; ++col) chain.accepted_at(0, col) = 1;

    std::vector<double> loglikprevious(nmax), loglikproposal(nmax), logpriorprevious(nmax), logpriorproposal(nmax);
    for (int n = 1; n <= nmax; ++n) {
        loglikprevious[n - 1] =
            models[n]->loglik(tauCprevious[n - 1], tauAprevious[n - 1], epsprevious[n - 1], ws[n - 1]);
    }
    for (double ll : loglikprevious) {
        if (ll == kNegInf) throw std::runtime_error("The data does not support the initial value.");
    }
    for (int n = 1; n <= nmax; ++n) {
        logpriorprevious[n - 1] = logprior(tauCprevious[n - 1], tauAprevious[n - 1], epsprevious[n - 1]);
    }
    double logpdfnprevious = s.prioronn.logpdf(nprevious);
    for (double lp : logpriorprevious) {
        if (lp == kNegInf) throw std::runtime_error("The initial value is outside the support of the prior.");
        if (!std::isfinite(lp)) throw std::runtime_error("The prior density is not finite at the initial value.");
    }
    if (logpdfnprevious == kNegInf) {
        throw std::runtime_error("The initial value is outside the support of the prior on n.");
    }
    for (int n = 1; n <= nmax; ++n) chain.logjointprob_at(0, n) = loglikprevious[n - 1] + logpriorprevious[n - 1];
    chain.logjointprob_at(0, 0) =
        loglikprevious[nprevious - 1] + logpriorprevious[nprevious - 1] + logpdfnprevious;

    // Proposal distributions: independent normal steps for tauC, tauA (variance 1e-6) and eps (sd 1e-4).
    std::vector<double> sdTauC(nmax, std::sqrt(1.0e-6)), sdTauA(nmax, std::sqrt(1.0e-6)), sdEps(nmax, 1.0e-4);
    std::vector<double> scalefactor(nmax, 2.38 * 2.38 / 4);  // 2.38^2/d with d = 4 parameters
    std::vector<double> tauCproposal(nmax), tauAproposal(nmax), epsproposal(nmax);
    const std::function<void(std::size_t)> evaluate_proposal = [&](std::size_t task) {
        const int i = nmax - 1 - static_cast<int>(task);
        loglikproposal[i] = models[i + 1]->loglik(tauCproposal[i], tauAproposal[i], epsproposal[i], ws[i]);
        logpriorproposal[i] = logprior(tauCproposal[i], tauAproposal[i], epsproposal[i]);
    };
    // With 80% probability the proposal for n is the current n, otherwise one of the other values.
    const double otherprob = nmax > 1 ? 0.2 / (nmax - 1) : 0.0;

    for (int step = 2; step <= nsteps; ++step) {  // 1-based step numbers as in Charon.jl
        const int row = step - 1;
        for (int i = 0; i < nmax; ++i) {
            tauCproposal[i] = tauCprevious[i] + sdTauC[i] * rng.normal();
            tauAproposal[i] = tauAprevious[i] + sdTauA[i] * rng.normal();
            epsproposal[i] = epsprevious[i] + sdEps[i] * rng.normal();
        }
        team.run(static_cast<std::size_t>(nmax), evaluate_proposal);
        for (int n = 1; n <= nmax; ++n) {
            const int i = n - 1;
            const double logacceptance =
                loglikproposal[i] - loglikprevious[i] + logpriorproposal[i] - logpriorprevious[i];
            if (rng.uniform() < std::exp(logacceptance)) {
                chain.tauC_at(row, n) = tauCprevious[i] = tauCproposal[i];
                chain.tauA_at(row, n) = tauAprevious[i] = tauAproposal[i];
                chain.eps_at(row, n) = epsprevious[i] = epsproposal[i];
                chain.accepted_at(row, n) = 1;
                chain.logjointprob_at(row, n) = loglikproposal[i] + logpriorproposal[i];
                loglikprevious[i] = loglikproposal[i];
                logpriorprevious[i] = logpriorproposal[i];
            } else {
                chain.tauC_at(row, n) = tauCprevious[i];
                chain.tauA_at(row, n) = tauAprevious[i];
                chain.eps_at(row, n) = epsprevious[i];
                chain.accepted_at(row, n) = 0;
                chain.logjointprob_at(row, n) = loglikprevious[i] + logpriorprevious[i];
            }
        }

        // The chain for n.
        int nproposal = nprevious;
        if (nmax > 1) {
            const double u = rng.uniform();
            double cumulative = 0.0;
            nproposal = nmax;
            for (int k = 1; k <= nmax; ++k) {
                cumulative += k == nprevious ? 0.8 : otherprob;
                if (u < cumulative) {
                    nproposal = k;
                    break;
                }
            }
        }
        if (nproposal == nprevious) {
            chain.nsample[row] = nproposal;
            chain.accepted_at(row, 0) = 1;
            chain.logjointprob_at(row, 0) = chain.logjointprob_at(row, nproposal) + logpdfnprevious;
        } else {
            const double logpdfnproposal = s.prioronn.logpdf(nproposal);
            const double logacceptprobn = chain.logjointprob_at(row, nproposal) -
                                          chain.logjointprob_at(row, nprevious) + logpdfnproposal - logpdfnprevious;
            if (rng.uniform() < std::exp(logacceptprobn)) {
                chain.nsample[row] = nprevious = nproposal;
                chain.accepted_at(row, 0) = 1;
                chain.logjointprob_at(row, 0) = chain.logjointprob_at(row, nproposal) + logpdfnproposal;
                logpdfnprevious = logpdfnproposal;
            } else {
                chain.nsample[row] = nprevious;
                chain.accepted_at(row, 0) = 0;
                chain.logjointprob_at(row, 0) = chain.logjointprob_at(row, nprevious) + logpdfnprevious;
            }
        }

        // Adaptation: proposal covariances from the thinned history (steps 3000:20:step).
        if (step % 1000 == 0 && step >= 5000) {
            for (int n = 1; n <= nmax; ++n) {
                const int i = n - 1;
                double mean = 0.0, varTauC = 0.0, varTauA = 0.0, varEps = 0.0;
                strided_moments(chain.tauC, nmax, i, 2999, step - 1, 20, mean, varTauC);
                strided_moments(chain.tauA, nmax, i, 2999, step - 1, 20, mean, varTauA);
                strided_moments(chain.eps, nmax, i, 2999, step - 1, 20, mean, varEps);
                if (10000 < step && step < nsteps / 4) {
                    int acceptedcount = 0;
                    for (int r = step - 1000; r <= step - 1; ++r) acceptedcount += chain.accepted_at(r, n);
                    const double acceptancerate = acceptedcount / 1000.0;  // last 1000 steps
                    if (acceptancerate < 0.1) {  // too low: take smaller steps
                        scalefactor[i] = 3 * scalefactor[i] / 4;
                        if (s.scalingmessages) {
                            log.info("Chain " + std::to_string(chainindex) +
                                     " talking, I scaled the scalefactor with factor 3/4, new scalefactor = " +
                                     format_vector(scalefactor) + ".");
                        }
                    } else if (acceptancerate > 0.15) {  // too high: take bigger steps
                        scalefactor[i] = 4 * scalefactor[i] / 3;
                        if (s.scalingmessages) {
                            log.info("Chain " + std::to_string(chainindex) +
                                     " talking, I scaled the scalefactor with factor 4/3, new scalefactor = " +
                                     format_vector(scalefactor) + ".");
                        }
                    }
                }
                sdTauC[i] = std::sqrt(scalefactor[i] * varTauC + 1.0e-14);
                sdTauA[i] = std::sqrt(scalefactor[i] * varTauA + 1.0e-14);
                sdEps[i] = std::sqrt(scalefactor[i]) * std::sqrt(varEps) + 1.0e-14;
            }
        }
        if (messages > 0 && step % messages == 0) {
            log.info("Chain " + std::to_string(chainindex) + " talking, I finished step " + std::to_string(step) +
                     " of " + std::to_string(nsteps) + ", " +
                     std::to_string(static_cast<int>(std::floor(100.0 * step / nsteps))) + " percent finished.");
        }
    }
}

}  // namespace

int resolve_threads(int requested) {
    if (requested > 0) return requested;
    const unsigned hw = std::thread::hardware_concurrency();
    return hw > 0 ? static_cast<int>(hw) : 1;
}

Chain::Chain(int nsteps_, int nmax_)
    : nsteps(nsteps_),
      nmax(nmax_),
      nsample(static_cast<std::size_t>(nsteps_)),
      tauC(static_cast<std::size_t>(nsteps_) * nmax_),
      tauA(static_cast<std::size_t>(nsteps_) * nmax_),
      eps(static_cast<std::size_t>(nsteps_) * nmax_),
      accepted(static_cast<std::size_t>(nsteps_) * (nmax_ + 1)),
      logjointprob(static_cast<std::size_t>(nsteps_) * (nmax_ + 1)) {}

GridResult exactposterior(const PreparedData& data, const std::vector<int>& nrange, const std::vector<double>& tauCrange,
                          const std::vector<double>& tauArange, const std::vector<double>& epsrange, int threads,
                          std::int64_t messages, Logger& log) {
    check_grid(nrange, tauCrange, tauArange, epsrange);
    threads = resolve_threads(threads);
    const ModelTable models = build_models(data, nrange, threads);
    return evaluate_grid(models, nrange, tauCrange, tauArange, epsrange, threads, messages, log);
}

std::vector<Chain> mcmc_sampler(const PreparedData& data, const MCMCSettings& s, Logger& log) {
    // The checks of Charon.jl.
    if (s.prioronn.minimum() < 1) {
        throw std::invalid_argument(
            "prioronn should have support on positive integers. There is always at least one individual.");
    }
    const int nmax = s.prioronn.maximum();
    if (nmax > kMaxIndividualsMCMC) {
        throw std::invalid_argument(
            "For the current implementation we require that the prior on n is bounded by 10. I.e. P(prioronn ≤ 10)=1");
    }
    if (s.prioronEps.lower() < 0.0) {
        throw std::invalid_argument("prioronϵ should have support on [0,0.5], but " + s.prioronEps.describe() +
                                    " can take negative values.");
    }
    if (s.prioronEps.upper() > 0.5) {
        throw std::invalid_argument("prioronϵ should have support on [0,0.5], but " + s.prioronEps.describe() +
                                    " can take values larger than 0.5.");
    }
    if (s.prioronTauC.lower() < 0.0 || s.prioronTauA.lower() < 0.0) {
        throw std::invalid_argument("The priors on τC and τA should have support on non-negative numbers.");
    }
    if (s.nchains <= 0) throw std::invalid_argument("There should be at least one chain.");
    if (s.nsteps <= 0) throw std::invalid_argument("nsteps should be positive.");

    const int threads = resolve_threads(s.threads);
    const std::int64_t messages = s.messages < 0 ? s.nsteps / 100 : s.messages;

    // Starting values: maximum likelihood on a grid that avoids the boundaries 0 (and 0.5 for ϵ).
    std::vector<int> nrange;
    for (int n = 1; n <= nmax; ++n) nrange.push_back(n);
    std::vector<double> gridTauC, gridTauA, gridEps;
    for (int i = 1; i <= 22; ++i) {
        gridTauC.push_back(s.prioronTauC.lower() + (s.prioronTauC.upper() - s.prioronTauC.lower()) * i / 22.0);
        gridTauA.push_back(s.prioronTauA.lower() + (s.prioronTauA.upper() - s.prioronTauA.lower()) * i / 22.0);
        gridEps.push_back(s.prioronEps.lower() + (s.prioronEps.upper() - s.prioronEps.lower()) * i / 23.0);
    }
    log.info("Building the model for n = 1, ..., " + std::to_string(nmax) + ".");
    const ModelTable models = build_models(data, nrange, threads);
    log.info("Searching starting values on a grid of " + std::to_string(nmax * 22 * 22 * 22) + " points.");
    const GridResult grid = evaluate_grid(models, nrange, gridTauC, gridTauA, gridEps, threads, 0, log);

    StartingValues start;
    start.tauC.resize(nmax);
    start.tauA.resize(nmax);
    start.eps.resize(nmax);
    std::vector<double> best(static_cast<std::size_t>(nmax) + 1, kNegInf);
    std::vector<bool> seen(static_cast<std::size_t>(nmax) + 1, false);
    double bestoverall = kNegInf;
    bool foundoverall = false;
    for (std::size_t i = 0; i < grid.loglik.size(); ++i) {  // first maximum, like Julia's argmax
        const int n = grid.n[i];
        const double ll = grid.loglik[i];
        // Skip grid points where the prior density is infinite (e.g. a beta prior with
        // B < 1 at its upper bound, which is on the grid): a chain could never leave them.
        const double lp = s.prioronTauC.logpdf(grid.tauC[i]) + s.prioronTauA.logpdf(grid.tauA[i]) +
                          s.prioronEps.logpdf(grid.eps[i]);
        if (!std::isfinite(lp)) continue;
        if (!seen[n] || ll > best[n]) {
            seen[n] = true;
            best[n] = ll;
            start.tauC[n - 1] = grid.tauC[i];
            start.tauA[n - 1] = grid.tauA[i];
            start.eps[n - 1] = grid.eps[i];
        }
        // The starting n maximises the likelihood among the values the prior allows.
        if (s.prioronn.logpdf(n) > kNegInf && (!foundoverall || ll > bestoverall)) {
            foundoverall = true;
            bestoverall = ll;
            start.n = n;
        }
    }
    for (int n = 1; n <= nmax; ++n) {
        if (!seen[n]) throw std::runtime_error("The prior density is not finite at any of the grid points for the starting values.");
    }
    log.info("Starting value for n: " + std::to_string(start.n) + ".");

    std::vector<Chain> chains;
    chains.reserve(s.nchains);
    for (int c = 0; c < s.nchains; ++c) chains.emplace_back(s.nsteps, nmax);
    std::vector<std::uint64_t> seeds(s.nchains);
    std::uint64_t state = s.seed;
    for (auto& seed : seeds) seed = splitmix64(state);
    // Chains run in parallel; spare threads help within the chains.
    const int concurrent = std::min(threads, s.nchains);
    const int perchain = std::max(1, threads / concurrent);
    log.info("Running " + std::to_string(s.nchains) + " chain(s) of " + std::to_string(s.nsteps) + " steps on " +
             std::to_string(concurrent * std::min(perchain, nmax)) + " thread(s).");
    parallel_for(static_cast<std::size_t>(s.nchains), concurrent, [&](std::size_t c) {
        run_chain(static_cast<int>(c) + 1, chains[c], models, start, s, messages, seeds[c], perchain, log);
    });
    return chains;
}

void write_unpacked_posterior(const std::vector<Chain>& chains, OutputFile& out) {
    out.write("nsample,τCsample,τAsample,ϵsample,accepted,logjointprob,chainid\n");
    std::string line;
    for (std::size_t c = 0; c < chains.size(); ++c) {
        const Chain& ch = chains[c];
        const std::string chainid = std::to_string(c + 1);
        double prevC = 0.0, prevA = 0.0, prevE = 0.0;
        int prevN = 0;
        for (int step = 0; step < ch.nsteps; ++step) {
            const int n = ch.nsample[step];
            const std::size_t idx = static_cast<std::size_t>(step) * ch.nmax + (n - 1);
            const double tC = ch.tauC[idx], tA = ch.tauA[idx], e = ch.eps[idx];
            // accepted: did any of n, τC, τA or ϵ change compared with the previous step?
            const bool accepted = step == 0 || !(tC == prevC && tA == prevA && e == prevE && n == prevN);
            line.clear();
            line += std::to_string(n);
            line += ',';
            line += format_number(tC);
            line += ',';
            line += format_number(tA);
            line += ',';
            line += format_number(e);
            line += accepted ? ",true," : ",false,";
            line += format_number(ch.logjointprob[static_cast<std::size_t>(step) * (ch.nmax + 1)]);
            line += ',';
            line += chainid;
            line += '\n';
            out.write(line);
            prevC = tC;
            prevA = tA;
            prevE = e;
            prevN = n;
        }
    }
}

void write_raw_chain(const Chain& ch, OutputFile& out) {
    std::string header = "nsample";
    for (const char* name : {"τC", "τA", "ϵ"}) {
        for (int n = 1; n <= ch.nmax; ++n) header += std::string(",") + name + "_" + std::to_string(n);
    }
    header += ",accepted_nsample";
    for (int n = 1; n <= ch.nmax; ++n) header += ",accepted_" + std::to_string(n);
    header += ",logjointprob";
    for (int n = 1; n <= ch.nmax; ++n) header += ",logjointprob_" + std::to_string(n);
    out.write(header + "\n");
    std::string line;
    for (int step = 0; step < ch.nsteps; ++step) {
        const std::size_t row = static_cast<std::size_t>(step) * ch.nmax;
        const std::size_t row1 = static_cast<std::size_t>(step) * (ch.nmax + 1);
        line = std::to_string(ch.nsample[step]);
        for (const std::vector<double>* m : {&ch.tauC, &ch.tauA, &ch.eps}) {
            for (int i = 0; i < ch.nmax; ++i) line += "," + format_number((*m)[row + i]);
        }
        for (int i = 0; i <= ch.nmax; ++i) line += ch.accepted[row1 + i] ? ",true" : ",false";
        for (int i = 0; i <= ch.nmax; ++i) line += "," + format_number(ch.logjointprob[row1 + i]);
        line += '\n';
        out.write(line);
    }
}

void write_grid(const GridResult& g, OutputFile& out) {
    out.write("n,τC,τA,ϵ,loglik\n");
    std::string line;
    for (std::size_t i = 0; i < g.loglik.size(); ++i) {
        line = std::to_string(g.n[i]) + "," + format_number(g.tauC[i]) + "," + format_number(g.tauA[i]) + "," +
               format_number(g.eps[i]) + "," + format_number(g.loglik[i]) + "\n";
        out.write(line);
    }
}

}  // namespace charon
