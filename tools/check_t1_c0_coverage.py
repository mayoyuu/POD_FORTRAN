#!/usr/bin/env python3
"""Check that the T1 ensemble covers the SRP-free C0 moments hour by hour.

For every hour of a campaign the script compares <orbit>/T1/<method> against
<orbit>/C0/<method> (hourly_moments.csv from the unified driver) and reports

  * m6, the 6D Mahalanobis distance of the C0 mean from the T1 mean measured in
    C0 covariance units,
  * lambda_min, the smallest eigenvalue of the C0-whitened difference
    Sigma_C0^-1/2 (Sigma_T1 - Sigma_C0) Sigma_C0^-1/2, which is dimensionless:
    0 means equal variance, +1 means T1 doubles the C0 variance along that
    mode and a negative value means T1 is tighter than C0 there,
  * the covariance volume ratio det(Sigma_T1)/det(Sigma_C0) and the smallest
    per-axis variance ratio.

A T1 prior that reaches SRP~0 must contain the C0 cloud: the mean shift stays
within a few C0 sigma and Sigma_T1 - Sigma_C0 stays positive semidefinite up to
Monte-Carlo noise. Exit status is 0 only when every hour, orbit and method
passes; the worst hour per metric is reported with the tolerances in use.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys
import numpy as np

MEANS = tuple(f'mean{i}' for i in range(1, 7))
COVARIANCE = tuple(f'cov{i}{j}' for i in range(1, 7) for j in range(1, 7))
VALID = ('accuracy_valid', 'precision_confirmed')


def sym(a):
    return (a+a.T)/2.


def spd_inverse(a, floor_rel):
    """Inverse of a symmetric matrix with eigenvalues floored at floor_rel*max."""
    w, v = np.linalg.eigh(sym(a))
    top = float(np.max(w))
    if not math.isfinite(top) or top <= 0.:
        raise ValueError('covariance has no positive eigenvalue')
    w = np.where(w > floor_rel*top, w, np.inf)
    return (v*(1./w)) @ v.T


def whitener(a, floor_rel):
    """Symmetric square root of the inverse, so W @ a @ W is the identity."""
    w, v = np.linalg.eigh(sym(a))
    top = float(np.max(w))
    if not math.isfinite(top) or top <= 0.:
        raise ValueError('covariance has no positive eigenvalue')
    w = np.where(w > floor_rel*top, w, floor_rel*top)
    return (v/np.sqrt(w)) @ v.T


def logdet(a):
    w = np.linalg.eigvalsh(sym(a))
    return float(np.sum(np.log(w))) if np.all(w > 0.) else float('-inf')


def read_moments(path):
    if not Path(path).is_file():
        raise FileNotFoundError(path)
    table = np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True, dtype=float))
    names = table.dtype.names
    missing = [k for k in ('hour', *MEANS, *COVARIANCE) if k not in names]
    if missing:
        raise ValueError(f'{path} is missing columns {missing}')
    hours = np.asarray(table['hour'], dtype=float)
    means = np.column_stack([table[k] for k in MEANS])
    cov = np.column_stack([table[k] for k in COVARIANCE]).reshape(-1, 6, 6)
    valid = np.ones(len(table), dtype=bool)
    for key in VALID:
        if key in names:
            valid &= np.asarray(table[key], dtype=float) == 1.
    return hours, means, cov, valid


def hour_metrics(mu0, s0, mu1, s1, floor_rel=1e-12):
    """Distance and covariance-containment metrics for one hour."""
    delta = mu1-mu0
    m6 = math.sqrt(max(float(delta @ spd_inverse(s0, floor_rel) @ delta), 0.))
    # Whitened difference: +1 means Sigma_T1 doubles the C0 variance along that
    # mode, 0 means equal, negative means T1 is tighter than C0 there.
    root = whitener(s0, floor_rel)
    difference = np.linalg.eigvalsh(sym(root @ sym(s1-s0) @ root))
    volume = logdet(s1)-logdet(s0)
    ratio = math.exp(min(volume, 700.)) if math.isfinite(volume) else float('inf')
    axes = np.diag(sym(s1))/np.diag(sym(s0))
    return dict(m6=m6, lambda_min=float(difference.min()),
                lambda_relative=float(difference.min()),
                volume_ratio=ratio, min_variance_ratio=float(np.min(axes)),
                mean_shift=float(math.sqrt(float(delta @ delta))))


def check_pair(root, orbit, method, args):
    base = root/orbit
    c0_path, t1_path = base/'C0'/method/'hourly_moments.csv', base/'T1'/method/'hourly_moments.csv'
    try:
        h0, m0, s0, v0 = read_moments(c0_path)
        h1, m1, s1, v1 = read_moments(t1_path)
    except (FileNotFoundError, ValueError) as exc:
        return dict(orbit=orbit, method=method, passed=False, reason=str(exc), hours=[])
    if len(h0) != len(h1) or not np.allclose(h0, h1, rtol=0, atol=1e-8):
        return dict(orbit=orbit, method=method, passed=False, hours=[],
                    reason=f'hour grids differ: C0 has {len(h0)} rows, T1 has {len(h1)}')
    keep = v0 & v1
    if np.count_nonzero(keep) < 2:
        return dict(orbit=orbit, method=method, passed=False, hours=[],
                    reason=f'only {np.count_nonzero(keep)} hours are numerically valid in both cases')
    hours = []
    for index in np.flatnonzero(keep):
        metrics = hour_metrics(m0[index], s0[index], m1[index], s1[index], args.relative_floor)
        hours.append(dict(hour=float(h0[index]), **metrics))
    worst_m6 = max(hours, key=lambda h: h['m6'])
    worst_lambda = min(hours, key=lambda h: h['lambda_relative'])
    worst_volume = min(hours, key=lambda h: h['volume_ratio'])
    worst_axis = min(hours, key=lambda h: h['min_variance_ratio'])
    checks = [
        dict(name='mean shift stays inside the C0 cloud', ok=bool(worst_m6['m6'] <= args.max_m6),
             detail=f"max m6={worst_m6['m6']:.3f} at h={worst_m6['hour']:g} (limit {args.max_m6:g})"),
        dict(name='Sigma_T1 - Sigma_C0 is positive semidefinite',
             ok=bool(worst_lambda['lambda_relative'] >= args.lambda_tol),
             detail=f"min lambda_rel={worst_lambda['lambda_relative']:+.5f} at h={worst_lambda['hour']:g} "
                    f'(limit {args.lambda_tol:+.3f})'),
        dict(name='T1 covariance volume is not smaller than C0',
             ok=bool(worst_volume['volume_ratio'] >= args.volume_tol),
             detail=f"min volume ratio={worst_volume['volume_ratio']:.4f} at h={worst_volume['hour']:g} "
                    f'(limit {args.volume_tol:g})'),
    ]
    return dict(orbit=orbit, method=method, passed=all(c['ok'] for c in checks), checks=checks,
                hour_count=len(hours), skipped_hours=int(len(h0)-len(hours)),
                worst=dict(m6=worst_m6, lambda_relative=worst_lambda, volume_ratio=worst_volume,
                           min_variance_ratio=worst_axis),
                hours=hours)


def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root', type=Path, required=True, help='campaign output root, e.g. ~/SRP_effect')
    p.add_argument('--orbits', nargs='+', help='orbit subdirectories (default: every orbit with C0 and T1)')
    p.add_argument('--methods', '--method', nargs='+', default=['UT'])
    p.add_argument('--max-m6', type=float, default=6.,
                   help='largest acceptable 6D Mahalanobis distance of the C0 mean (default 6)')
    p.add_argument('--lambda-tol', type=float, default=-.05,
                   help='allowed negative whitened eigenvalue of Sigma_T1-Sigma_C0 (-0.05 = 5%% tighter than C0)')
    p.add_argument('--volume-tol', type=float, default=.9,
                   help='smallest acceptable det(Sigma_T1)/det(Sigma_C0) (default 0.9)')
    p.add_argument('--relative-floor', type=float, default=1e-12,
                   help='relative eigenvalue floor used when inverting Sigma_C0 (default 1e-12)')
    p.add_argument('--json', type=Path, help='write the full per-hour report here')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'campaign root does not exist: {root}')
    orbits = args.orbits or sorted({p.parent.parent.parent.name
                                    for p in root.glob('*/T1/*/hourly_moments.csv')})
    if not orbits:
        raise ValueError(f'no <orbit>/T1/.../hourly_moments.csv under {root}')
    reports = [check_pair(root, orbit, method, args) for orbit in orbits for method in args.methods]
    for report in reports:
        if report.get('reason'):
            print(f"[{report['orbit']}/{report['method']}] FAIL {report['reason']}")
            continue
        worst = report['worst']
        print(f"[{report['orbit']}/{report['method']}] {'PASS' if report['passed'] else 'FAIL'} "
              f"h={report['hour_count']} (skipped {report['skipped_hours']}) "
              f"max m6={worst['m6']['m6']:.3f}@{worst['m6']['hour']:g} "
              f"min lambda_rel={worst['lambda_relative']['lambda_relative']:+.4f}"
              f"@{worst['lambda_relative']['hour']:g} "
              f"min vol={worst['volume_ratio']['volume_ratio']:.3f}@{worst['volume_ratio']['hour']:g}")
        for item in report['checks']:
            print(f"  {'PASS' if item['ok'] else 'FAIL'}  {item['name']}: {item['detail']}")
    passed = sum(1 for r in reports if r['passed'])
    print(f'{passed}/{len(reports)} orbit/method pairs pass the T1-vs-C0 coverage check')
    if args.json:
        args.json.write_text(json.dumps(dict(schema='t1_c0_coverage_v1', root=str(root),
                                             max_m6=args.max_m6, lambda_tol=args.lambda_tol,
                                             volume_tol=args.volume_tol, pairs=reports),
                                       indent=2)+'\n', encoding='utf-8')
    return 0 if passed == len(reports) else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
