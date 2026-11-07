#pragma once

// The likelihood of CHARON (port of the likelihood part of Charon.jl).
//
// Model (van Waaij et al., G3 2026; Schraiber 2018): with n diploid individuals,
// drift times tauC (anchor) and tauA (sample population), and error rate eps,
// a locus with coverage R, d derived reads and anchor frequency y has probability
//
//     sum_k q_k (e^{tauA Q} e^{tauC Qdown} h(y))_k,
//     q_k = binom(2n, k) * Binomial(R, p_k).pdf(d),  p_k = k(1-eps)/(2n) + (2n-k)eps/(2n),
//     h(y)_k = y^k (1-y)^(2n-k),                      k = 0, ..., 2n.

#include <cstdint>
#include <vector>

#include <Eigen/Dense>

namespace charon {

using Mat = Eigen::MatrixXd;
using RowMat = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using Vec = Eigen::VectorXd;

/// Decomposition M = P * diag(D) * Pinv, used for fast matrix exponentials:
/// e^{sM} = P * diag(e^{sD}) * Pinv.
struct EigenExpansion {
    Mat P;
    Vec D;
    Mat Pinv;

    EigenExpansion() = default;
    /// Throws std::runtime_error when M has complex eigenvalues.
    explicit EigenExpansion(const Mat& M);

    /// e^{sM}.
    Mat exp(double s) const;
    /// e^{sM} written into out; tmp is scratch space.
    void exp(double s, Mat& out, Mat& tmp) const;
};

/// The (2n+1)x(2n+1) tridiagonal matrix Q of Schraiber (2018).
Mat makeQ(int n);
/// The (2n+1)x(2n+1) tridiagonal matrix Qꜜ of Schraiber (2018).
Mat makeQdown(int n);
/// e^{tauA Q} e^{tauC Qdown}.
Mat calcmatrix(double tauC, double tauA, const EigenExpansion& Qdown, const EigenExpansion& Q);

/// The binomial coefficients N over k, k = 0, ..., N (exact for N <= 60).
std::vector<double> binomial_coefficients(int N);
/// The binomial probability mass function, P(Binomial(R, p) = d).
double binomial_pdf(std::int64_t R, double p, std::int64_t d);

/// The data after filtering and aggregation, shared by all n.
///
/// There are counts[i] loci with coverage pair_coverage[pair[i]], pair_derived[pair[i]]
/// derived reads and anchor frequency frequencies[i]. All coverages are positive and
/// all frequencies lie strictly between 0 and 1.
struct PreparedData {
    std::vector<std::int64_t> pair_coverage;  ///< coverage R of each distinct (R, d) pair
    std::vector<std::int64_t> pair_derived;   ///< derived reads d of each distinct (R, d) pair
    std::vector<std::int32_t> pair;           ///< per locus: index into the pair arrays
    std::vector<double> frequencies;          ///< per locus: anchor frequency y
    std::vector<double> counts;               ///< per locus: number of loci with this data

    std::size_t num_loci() const { return pair.size(); }
    std::size_t num_pairs() const { return pair_coverage.size(); }
    double total_count() const;
};

/// Everything needed to evaluate the likelihood conditioned on n individuals.
/// Immutable after construction, so one instance can be shared between threads;
/// every thread needs its own Workspace.
class ConditionalModel {
public:
    struct Workspace {
        RowMat q;  ///< num_pairs x (2n+1): q_k for every (R, d) pair
        RowMat B;  ///< num_pairs x (2n+1): q^T M for every (R, d) pair
        Mat expA, expC, M, tmp;
    };

    ConditionalModel(int n, const PreparedData& data);

    int n() const { return n_; }
    int dim() const { return m_; }
    Workspace make_workspace() const;

    /// The log-likelihood (logprobderivedreads! in Charon.jl). Returns -Inf when
    /// tauC < 0, tauA < 0, eps < 0 or eps > 1.
    double loglik(double tauC, double tauA, double eps, Workspace& ws) const;

    /// e^{tau Q} and e^{tau Qdown}, the building blocks of the matrix M.
    void exp_Q(double tau, Mat& out, Mat& tmp) const { Q_.exp(tau, out, tmp); }
    void exp_Qdown(double tau, Mat& out, Mat& tmp) const { Qdown_.exp(tau, out, tmp); }

    /// The log-likelihood for a given M = e^{tauA Q} e^{tauC Qdown}.
    double loglik_given_matrix(const Mat& M, double eps, Workspace& ws) const;

    /// Fill ws.q: q(pair, k) = binom(2n, k) P(Binomial(R, p_k) = d) (updateq! in Charon.jl).
    void update_q(double eps, Workspace& ws) const;

private:
    int n_;
    int m_;  // 2n + 1
    const PreparedData* data_;
    EigenExpansion Q_;
    EigenExpansion Qdown_;
    std::vector<double> binomial_;  // binom(2n, k)
    std::vector<double> lchoose_;   // log(binom(R, d)) per (R, d) pair
    std::vector<double> h_;         // num_loci x (2n+1), row-major: y^k (1-y)^(2n-k)
};

}  // namespace charon
