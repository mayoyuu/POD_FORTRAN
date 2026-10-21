"""Shared input design, augmented UT and Henze-Zirkler diagnostics.

Arrays are samples x variables. No scipy dependency; HZ uses the published
lognormal approximation and evaluates pair distances in bounded blocks.
"""
from __future__ import annotations
import math
import numpy as np


def shared_design(n: int, seed: int, sigma_limit: float = 3.0) -> np.ndarray:
    if n < 2 or not math.isfinite(sigma_limit) or sigma_limit <= 0:
        raise ValueError('sample count must be >=2 and sigma limit positive')
    state_rng, mode_rng = [np.random.Generator(np.random.PCG64DXSM(s))
                           for s in np.random.SeedSequence(seed).spawn(2)]
    chunks, left = [], n
    while left:
        z = state_rng.standard_normal((max(128, math.ceil(left*1.1)), 7))
        z = z[np.max(np.abs(z), axis=1) <= sigma_limit][:left]
        chunks.append(z)
        left -= len(z)
    return np.column_stack((np.concatenate(chunks), mode_rng.random((n, 2))))


def weighted_moments(points, wm, wc):
    points = np.asarray(points, dtype=float)
    wm, wc = np.asarray(wm, dtype=float), np.asarray(wc, dtype=float)
    if points.ndim != 2 or wm.shape != (len(points),) or wc.shape != wm.shape:
        raise ValueError('weight and point shapes do not match')
    if not np.all(np.isfinite(points)) or not np.all(np.isfinite([wm, wc])):
        raise ValueError('moments require finite inputs')
    if not np.isclose(wm.sum(), 1., rtol=1e-12, atol=1e-12):
        raise ValueError('mean weights must sum to one')
    # Anchor the sum near the cloud rather than subtracting huge orbital states.
    mean = points[0] + wm @ (points - points[0])
    centered = points - mean
    return mean, (centered.T * wc) @ centered


def joint_sigma_points(samples, alpha=.5, beta=2., kappa=0.):
    x = np.asarray(samples, dtype=float)
    if x.ndim != 2 or len(x) < 2 or x.shape[1] < 1 or not np.all(np.isfinite(x)):
        raise ValueError('UT needs a finite sample matrix with at least two rows')
    n = x.shape[1]
    c = alpha*alpha*(n+kappa)
    if not all(math.isfinite(v) for v in (alpha, beta, kappa)) or c <= 0:
        raise ValueError('invalid UT scaling parameters')
    mean = x.mean(axis=0)
    centered = x-mean
    std = np.sqrt(np.sum(centered*centered, axis=0)/(len(x)-1))
    scale = np.where(std > 0, std, 1.)
    normalized = centered/scale
    corr = normalized.T @ normalized/(len(x)-1)
    eig, vectors = np.linalg.eigh((corr+corr.T)/2)
    tol = 100*n*np.finfo(float).eps*max(1., np.max(eig))
    if np.min(eig) < -tol:
        raise ValueError('UT input covariance is indefinite')
    root = scale[:, None] * (vectors*np.sqrt(np.maximum(eig, 0.)))
    offsets = math.sqrt(c)*root.T
    points = np.vstack((mean, mean+offsets, mean-offsets))
    wm = np.full(2*n+1, 1./(2*c))
    wm[0] = 1-n/c
    wc = wm.copy()
    wc[0] += 1-alpha*alpha+beta
    return points, wm, wc


def henze_zirkler(samples, alpha=.05, block=256):
    x = np.asarray(samples, dtype=float)
    if x.ndim != 2 or len(x) <= x.shape[1]+1 or not np.all(np.isfinite(x)):
        raise ValueError('HZ needs finite observations and N > dimension+1')
    if not 0 < alpha < 1 or block < 1:
        raise ValueError('invalid HZ alpha/block')
    n, p = x.shape
    base = dict(n=n, dimension=p, alpha=alpha, statistic=None,
                p_value=None, reject_normal=None, rank=0)
    centered = x-x.mean(axis=0)
    std = np.std(centered, axis=0)
    if np.any(std == 0):
        return dict(base, status='singular_covariance')
    normalized = centered/std
    cov = normalized.T @ normalized/n
    w, v = np.linalg.eigh((cov+cov.T)/2)
    tol = max(w.max(), 1.)*1e-12
    rank = int(np.count_nonzero(w > tol))
    base['rank'] = rank
    if rank != p:
        return dict(base, status='singular_covariance')
    y = normalized @ (v/np.sqrt(w))
    di = np.einsum('ij,ij->i', y, y)
    b = ((2*p+1)*n/4.)**(1./(p+4))/math.sqrt(2.)
    b2, b4, b8 = b*b, b**4, b**8
    pairs = 0.
    for start in range(0, n, block):
        dist = di[start:start+block, None]+di[None, :]-2*y[start:start+block]@y.T
        pairs += np.exp(-b2/2*np.maximum(dist, 0.)).sum()
    hz = n*(pairs/n**2 - 2*(1+b2)**(-p/2)*
            np.exp(-b2/(2*(1+b2))*di).mean()+(1+2*b2)**(-p/2))
    a = 1+2*b2
    wb = (1+b2)*(1+3*b2)
    mu = 1-a**(-p/2)*(1+p*b2/a+p*(p+2)*b4/(2*a*a))
    variance = (2*(1+4*b2)**(-p/2)
                + 2*a**(-p)*(1+2*p*b4/a**2+3*p*(p+2)*b8/(4*a**4))
                - 4*wb**(-p/2)*(1+3*p*b4/(2*wb)+p*(p+2)*b8/(2*wb**2)))
    if hz <= 0 or mu <= 0 or variance <= 0:
        return dict(base, status='approximation_unavailable', statistic=float(hz))
    sigma_log = math.sqrt(math.log1p(variance/(mu*mu)))
    mean_log = math.log(mu)-sigma_log*sigma_log/2
    pvalue = .5*math.erfc((math.log(hz)-mean_log)/(math.sqrt(2)*sigma_log))
    return dict(base, status='ok', statistic=float(hz), p_value=pvalue,
                reject_normal=bool(pvalue < alpha))
