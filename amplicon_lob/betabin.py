"""Beta-binomial background model for the positional Limit of Blank.

The Gaussian LoB (mean + 1.645*SD) treats the per-sample blank fraction as normal. But the
blank is a *count* -- k non-reference reads of n -- and the between-sample spread is usually
larger than binomial sampling alone (batch/run effects): overdispersion. The beta-binomial
models this explicitly with a Beta(alpha, beta) on the per-site rate (the deepSNV/shearwater
family). We fit it by method of moments and take the 95th percentile of the fitted rate as the
noise ceiling.

Because the assay runs at very high depth, the finite-depth binomial sampling term is negligible
next to the rate spread, so the 95th percentile of the Beta rate distribution is an accurate LoB
for a newly observed VAF. The fit needs only three moments, so it can be computed either from the
raw per-sample counts (``betabin_lob``) or from stored per-site summary statistics
(``betabin_lob_from_moments``) -- both call the same estimator.
"""
import numpy as np
from scipy import stats


def betabin_lob_from_moments(m, var_p, inv_n_mean, q=0.95):
    """Beta-binomial LoB from summary moments of the per-sample blank rate.

    Parameters
    ----------
    m          : depth-weighted (or plain) mean blank rate across samples.
    var_p      : between-sample variance of the per-sample rates p_i = k_i / n_i.
    inv_n_mean : mean of 1 / n_i across samples (approx 1 / mean_depth when depths are similar).

    Returns ``(lob_vaf, rho)`` where rho is the intra-cluster correlation (overdispersion):
    0 = pure binomial sampling, ->1 = strong between-sample/batch spread.
    """
    if not np.isfinite(m) or m <= 0:
        return 0.0, 0.0
    if m >= 1:
        return 1.0, np.nan
    binom_var = m * (1.0 - m) * inv_n_mean          # variance expected from binomial sampling
    rho = (var_p - binom_var) / (m * (1.0 - m))     # excess spread -> intra-cluster correlation
    rho = float(min(max(rho, 1e-6), 0.999))         # clip to a valid (0, 1) Beta
    phi = (1.0 - rho) / rho                          # precision alpha + beta
    a, b = m * phi, (1.0 - m) * phi
    return float(stats.beta.ppf(q, a, b)), rho


def betabin_lob(k, n, q=0.95):
    """Beta-binomial LoB from per-sample counts (k alt of n reads). See module docstring."""
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)
    tot_n = n.sum()
    if len(k) < 2 or tot_n <= 0:
        return (k.sum() / tot_n if tot_n else 0.0), np.nan
    m = k.sum() / tot_n                             # depth-weighted pooled mean rate
    p = k / n
    return betabin_lob_from_moments(m, p.var(ddof=1), float(np.mean(1.0 / n)), q)
