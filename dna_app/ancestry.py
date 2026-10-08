"""Supervised admixture estimate against the 1000 Genomes superpopulations.

Allele frequencies for the five superpopulations are held fixed and only the
individual's mixture proportions ``q`` are estimated, with the same EM update
ADMIXTURE uses in supervised mode:

    p_mix(site)  = sum_k q_k * p_k(site)
    a_k(site)    = q_k * p_k / p_mix                 (alt allele responsibility)
    b_k(site)    = q_k * (1 - p_k) / (1 - p_mix)     (ref allele responsibility)
    q_k'         = sum_sites [g * a_k + (2 - g) * b_k] / (2 * n_sites)
"""

from __future__ import annotations

import numpy as np


def estimate_admixture(genotypes: np.ndarray, freqs: np.ndarray, iters: int = 500,
                       tol: float = 1e-7) -> tuple[np.ndarray, float]:
    """Return (q, log_likelihood).

    genotypes: shape (n,) alt-allele counts in {0, 1, 2}
    freqs:     shape (n, K) alt-allele frequency per population
    """
    g = genotypes.astype(np.float64)[:, None]
    p = np.clip(freqs.astype(np.float64), 1e-3, 1 - 1e-3)
    n, k = p.shape
    q = np.full(k, 1.0 / k)
    for _ in range(iters):
        mix = p @ q
        a = (p * q) / mix[:, None]
        b = ((1 - p) * q) / (1 - mix)[:, None]
        new_q = (g * a + (2 - g) * b).sum(axis=0) / (2 * n)
        new_q /= new_q.sum()
        if np.abs(new_q - q).max() < tol:
            q = new_q
            break
        q = new_q
    mix = np.clip(p @ q, 1e-12, 1 - 1e-12)
    ll = float((g[:, 0] * np.log(mix) + (2 - g[:, 0]) * np.log(1 - mix)).sum())
    return q, ll


def bootstrap_ci(genotypes: np.ndarray, freqs: np.ndarray, reps: int = 30,
                 seed: int = 0) -> np.ndarray:
    """95% bootstrap interval per population, shape (2, K)."""
    rng = np.random.default_rng(seed)
    n = len(genotypes)
    qs = []
    for _ in range(reps):
        idx = rng.integers(0, n, n)
        q, _ = estimate_admixture(genotypes[idx], freqs[idx], iters=200, tol=1e-5)
        qs.append(q)
    return np.percentile(np.array(qs), [2.5, 97.5], axis=0)
