#include "charon/model.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#include <Eigen/Eigenvalues>

namespace charon {

namespace {

constexpr double kNegInf = -std::numeric_limits<double>::infinity();

double log_binomial_coefficient(std::int64_t R, std::int64_t d) {
    return std::lgamma(static_cast<double>(R) + 1.0) - std::lgamma(static_cast<double>(d) + 1.0) -
           std::lgamma(static_cast<double>(R - d) + 1.0);
}

// P(Binomial(R, p) = d) given log(binom(R, d)), log(p) and log(1-p).
inline double binomial_pdf_from_logs(std::int64_t R, std::int64_t d, double p, double lchoose,
                                     double logp, double log1mp) {
    if (p <= 0.0) return d == 0 ? 1.0 : 0.0;
    if (p >= 1.0) return d == R ? 1.0 : 0.0;
    return std::exp(lchoose + static_cast<double>(d) * logp + static_cast<double>(R - d) * log1mp);
}

// The success probability p_k = k(1-eps)/(2n) + (2n-k)eps/(2n), as written in Charon.jl.
inline double success_probability(int k, int twon, double eps) {
    return k * (1.0 - eps) / twon + (twon - k) * eps / twon;
}

/// Balance M as LAPACK's dgebal does (the scaling step, as used by Julia's eigen):
/// returns powers of two d such that diag(d)^{-1} M diag(d) has rows and columns of
/// similar norm. The scaling is exact and makes the eigenvectors better conditioned.
Vec balancing_scale(Mat& M) {
    constexpr double radix = 2.0;
    const Eigen::Index m = M.rows();
    Vec d = Vec::Ones(m);
    bool converged = false;
    while (!converged) {
        converged = true;
        for (Eigen::Index i = 0; i < m; ++i) {
            double c = M.col(i).cwiseAbs().sum() - std::fabs(M(i, i));
            double r = M.row(i).cwiseAbs().sum() - std::fabs(M(i, i));
            if (c == 0.0 || r == 0.0) continue;
            // EISPACK balanc: find the power of two f that best equalises f*c and r/f.
            const double s = c + r;
            double f = 1.0;
            double g = r / radix;
            while (c < g) {
                f *= radix;
                c *= radix * radix;
            }
            g = r * radix;
            while (c >= g) {
                f /= radix;
                c /= radix * radix;
            }
            if ((c + r) / f < 0.95 * s) {
                converged = false;
                d(i) *= f;
                M.row(i) /= f;
                M.col(i) *= f;
            }
        }
    }
    return d;
}

}  // namespace

EigenExpansion::EigenExpansion(const Mat& input) {
    // Decompose the balanced matrix B = diag(d)^{-1} M diag(d): if B = Pb D Pb^{-1}, then
    // M = (diag(d) Pb) D (Pb^{-1} diag(d)^{-1}).
    Mat M = input;
    const Vec d = balancing_scale(M);
    Eigen::EigenSolver<Mat> solver(M, /*computeEigenvectors=*/true);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("EigenExpansion: the eigendecomposition did not converge.");
    }
    const Eigen::VectorXcd values = solver.eigenvalues();
    const Eigen::MatrixXcd vectors = solver.eigenvectors();
    const double scale = std::max(1.0, values.cwiseAbs().maxCoeff());
    if (values.imag().cwiseAbs().maxCoeff() > 1e-12 * scale ||
        vectors.imag().cwiseAbs().maxCoeff() > 1e-12) {
        throw std::runtime_error("EigenExpansion: the matrix has complex eigenvalues.");
    }
    const Mat Pb = vectors.real();
    D = values.real();
    P = d.asDiagonal() * Pb;
    Pinv = Pb.inverse() * d.cwiseInverse().asDiagonal();
}

Mat EigenExpansion::exp(double s) const {
    Mat out, tmp;
    exp(s, out, tmp);
    return out;
}

void EigenExpansion::exp(double s, Mat& out, Mat& tmp) const {
    // (P * e^{sD}) * Pinv, the same association as Charon.jl.
    tmp.resize(P.rows(), P.cols());
    for (Eigen::Index j = 0; j < P.cols(); ++j) tmp.col(j) = P.col(j) * std::exp(s * D(j));
    out.noalias() = tmp * Pinv;
}

Mat makeQ(int n) {
    const int m = 2 * n + 1;
    const double twon = 2.0 * n;
    Mat Q = Mat::Zero(m, m);
    for (int k = 0; k < m; ++k) {
        Q(k, k) = -k * (twon - k);
        if (k >= 1) Q(k, k - 1) = 0.5 * k * (k - 1.0);
        if (k + 1 < m) Q(k, k + 1) = 0.5 * (twon - k) * (twon - k - 1.0);
    }
    return Q;
}

Mat makeQdown(int n) {
    const int m = 2 * n + 1;
    const double twon = 2.0 * n;
    Mat Q = Mat::Zero(m, m);
    for (int k = 0; k < m; ++k) {
        Q(k, k) = -k * (twon - k + 1.0);
        if (k >= 1) Q(k, k - 1) = 0.5 * k * (k - 1.0);
        if (k + 1 < m) Q(k, k + 1) = 0.5 * (twon - k + 1.0) * (twon - k);
    }
    return Q;
}

Mat calcmatrix(double tauC, double tauA, const EigenExpansion& Qdown, const EigenExpansion& Q) {
    return Q.exp(tauA) * Qdown.exp(tauC);
}

std::vector<double> binomial_coefficients(int N) {
    std::vector<double> c(static_cast<std::size_t>(N) + 1);
    std::uint64_t value = 1;  // exact for N <= 60
    for (int k = 0; k <= N; ++k) {
        c[k] = static_cast<double>(value);
        value = value * static_cast<std::uint64_t>(N - k) / static_cast<std::uint64_t>(k + 1);
    }
    return c;
}

double binomial_pdf(std::int64_t R, double p, std::int64_t d) {
    if (d < 0 || d > R) return 0.0;
    return binomial_pdf_from_logs(R, d, p, log_binomial_coefficient(R, d), std::log(p), std::log1p(-p));
}

double PreparedData::total_count() const {
    double total = 0.0;
    for (double c : counts) total += c;
    return total;
}

ConditionalModel::ConditionalModel(int n, const PreparedData& data)
    : n_(n),
      m_(2 * n + 1),
      data_(&data),
      Q_(makeQ(n)),
      Qdown_(makeQdown(n)),
      binomial_(binomial_coefficients(2 * n)) {
    if (n < 1) throw std::invalid_argument("ConditionalModel: n should be positive.");
    // log(binom(R, d)) from a table of log factorials accumulated in long double. This
    // avoids std::lgamma, which is not thread-safe in glibc (it writes signgam).
    std::int64_t maxR = 0;
    for (std::int64_t R : data.pair_coverage) maxR = std::max(maxR, R);
    std::vector<long double> logfactorial(static_cast<std::size_t>(maxR) + 1, 0.0L);
    for (std::int64_t k = 2; k <= maxR; ++k) {
        logfactorial[k] = logfactorial[k - 1] + std::log(static_cast<long double>(k));
    }
    lchoose_.resize(data.num_pairs());
    for (std::size_t p = 0; p < data.num_pairs(); ++p) {
        const std::int64_t R = data.pair_coverage[p], d = data.pair_derived[p];
        lchoose_[p] = static_cast<double>(logfactorial[R] - logfactorial[d] - logfactorial[R - d]);
    }
    const std::size_t L = data.num_loci();
    h_.resize(L * static_cast<std::size_t>(m_));
    for (std::size_t i = 0; i < L; ++i) {
        const double y = data.frequencies[i];
        for (int k = 0; k < m_; ++k) {
            h_[i * m_ + k] = std::pow(y, k) * std::pow(1.0 - y, 2 * n - k);
        }
    }
}

ConditionalModel::Workspace ConditionalModel::make_workspace() const {
    Workspace ws;
    const auto P = static_cast<Eigen::Index>(data_->num_pairs());
    ws.q.resize(P, m_);
    ws.B.resize(P, m_);
    ws.expA.resize(m_, m_);
    ws.expC.resize(m_, m_);
    ws.M.resize(m_, m_);
    ws.tmp.resize(m_, m_);
    return ws;
}

void ConditionalModel::update_q(double eps, Workspace& ws) const {
    const int twon = 2 * n_;
    const std::size_t P = data_->num_pairs();
    for (int k = 0; k < m_; ++k) {
        const double p = success_probability(k, twon, eps);
        const double logp = std::log(p);
        const double log1mp = std::log1p(-p);
        for (std::size_t j = 0; j < P; ++j) {
            ws.q(static_cast<Eigen::Index>(j), k) =
                binomial_[k] * binomial_pdf_from_logs(data_->pair_coverage[j], data_->pair_derived[j], p,
                                                      lchoose_[j], logp, log1mp);
        }
    }
}

double ConditionalModel::loglik_given_matrix(const Mat& M, double eps, Workspace& ws) const {
    update_q(eps, ws);
    // B = q^T M for every (R, d) pair, so that each locus only needs a dot product with h(y).
    // This is the association order of Julia's dot(q, M, h).
    ws.B.noalias() = ws.q * M;
    const std::size_t L = data_->num_loci();
    const std::int32_t* pair = data_->pair.data();
    const double* counts = data_->counts.data();
    const double* h = h_.data();
    const double* B = ws.B.data();
    const int m = m_;
    double loglik = 0.0;
    for (std::size_t i = 0; i < L; ++i) {
        const double* b = B + static_cast<std::size_t>(pair[i]) * m;
        const double* hi = h + i * m;
        double s = 0.0;
        for (int k = 0; k < m; ++k) s += b[k] * hi[k];
        loglik += counts[i] * std::log(std::max(0.0, s));
    }
    return loglik;
}

double ConditionalModel::loglik(double tauC, double tauA, double eps, Workspace& ws) const {
    if (tauC < 0.0 || tauA < 0.0 || eps < 0.0 || eps > 1.0) return kNegInf;
    Q_.exp(tauA, ws.expA, ws.tmp);
    Qdown_.exp(tauC, ws.expC, ws.tmp);
    ws.M.noalias() = ws.expA * ws.expC;
    return loglik_given_matrix(ws.M, eps, ws);
}

}  // namespace charon
