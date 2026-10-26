#!/usr/bin/env python3
"""Audit the T1 eta prior actually written into a campaign directory.

Reads <root>/<orbit>/<case>/population.json and <case>/initial_points.csv and
verifies, without rerunning the propagation:

  * the recorded prior name/support matches what was requested,
  * every T1 eta lies inside the physical SRP bound [-1, 1] and the ADS domain,
  * for --expect-prior uniform, eta is statistically uniform on the domain
    (Kolmogorov-Smirnov against U(-a, +a) plus a std consistency band) and the
    SRP~0 tail is populated at the rate the uniform prior implies. Because the
    SRP scale factor is 1+eta, SRP~0 sits at eta=-1, not eta=0,
  * C0..C3 keep eta exactly zero and share one bit-identical 6D state design
    with T1 (input_state_sha256 equality), so the prior switch cannot have
    changed anything except the T1 eta column.

Exit status is 0 only when every check passes.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
import sys
import numpy as np

CASES = ('C0', 'C1', 'C2', 'C3', 'T1')


def load_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def read_points(path):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    return np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True, dtype=float))


def kolmogorov_sf(lam):
    """Asymptotic survival function of the Kolmogorov distribution."""
    if lam < 0.2:
        return 1.
    total = 0.
    for k in range(1, 201):
        term = (-1)**(k-1)*math.exp(-2.*k*k*lam*lam)
        total += term
        if abs(term) < 1e-13:
            break
    return min(1., max(0., 2.*total))


def uniform_ks(values, a):
    """Two-sided KS statistic of values against U(-a, a) and its p-value."""
    x = np.sort(np.asarray(values, dtype=float))
    n = len(x)
    cdf = np.clip((x+a)/(2.*a), 0., 1.)
    i = np.arange(1, n+1)
    d = max(float(np.max(i/n-cdf)), float(np.max(cdf-(i-1)/n)))
    lam = (math.sqrt(n)+0.12+0.11/math.sqrt(n))*d
    return d, kolmogorov_sf(lam)


def check(checks, name, ok, detail=''):
    checks.append(dict(name=name, ok=bool(ok), detail=detail))
    return bool(ok)


def check_orbit(root, orbit, expect_prior, eta_sigma, domain_sigma, alpha):
    checks, stats = [], dict(orbit=orbit)
    cases = {}
    for name in CASES:
        pop = root/orbit/name/'population.json'
        if not pop.is_file():
            check(checks, f'{name}/population.json', False, f'missing {pop}')
            continue
        cases[name] = load_json(pop)
    if 'T1' not in cases:
        return dict(orbit=orbit, passed=False, checks=checks, statistics=stats)

    t1 = cases['T1']
    recorded_support = t1.get('eta_support')
    support = float(recorded_support) if recorded_support is not None else domain_sigma*eta_sigma
    stats.update(eta_prior=t1.get('eta_prior'), eta_support=recorded_support,
                 eta_sigma=eta_sigma, domain_sigma=domain_sigma, support=float(support))
    check(checks, 'T1 records eta_prior', t1.get('eta_prior') is not None,
          f"eta_prior={t1.get('eta_prior')!r}")
    check(checks, 'T1 prior name matches request', t1.get('eta_prior') == expect_prior,
          f"recorded={t1.get('eta_prior')!r} expected={expect_prior!r}")
    check(checks, 'T1 records a positive eta_support',
          recorded_support is not None and support > 0., f'recorded={recorded_support}')
    if eta_sigma is not None:
        check(checks, 'T1 eta_support matches domain_sigma*eta_sigma',
              abs(support-domain_sigma*eta_sigma) <= 1e-12*max(1., support),
              f'recorded={support} expected={domain_sigma*eta_sigma}')
    check(checks, 'ADS domain stays inside the physical SRP bound', support <= 1.+1e-12,
          f'support={support}')

    try:
        t1_points = read_points(root/orbit/'T1'/'initial_points.csv')
    except FileNotFoundError as exc:
        check(checks, 'T1/initial_points.csv', False, str(exc))
        return dict(orbit=orbit, passed=False, checks=checks, statistics=stats)

    eta = np.asarray(t1_points['eta'], dtype=float)
    check(checks, 'T1 eta is finite', np.all(np.isfinite(eta)), f'n={len(eta)}')
    check(checks, 'T1 eta inside [-1, 1]',
          len(eta) > 0 and np.min(eta) >= -1. and np.max(eta) <= 1.,
          f'min={np.min(eta):.6e} max={np.max(eta):.6e}')
    if support > 0:
        check(checks, 'T1 eta inside the ADS domain',
              np.max(np.abs(eta)) <= support*(1.+1e-12),
              f'max|eta|={np.max(np.abs(eta)):.6e} support={support:.6e}')
    n = len(eta)
    mean = float(np.mean(eta))
    std = float(np.std(eta, ddof=1)) if n > 1 else float('nan')
    # The SRP scale factor is 1+eta, so SRP~0 sits at eta=-1, not eta=0.
    near_zero = float(np.mean(1.+eta <= .01))
    stats.update(sample_count=n, eta_min=float(np.min(eta)), eta_max=float(np.max(eta)),
                 eta_mean=mean, eta_std=std,
                 fraction_srp_scale_le_1pct=near_zero,
                 fraction_abs_eta_le_1pct=float(np.mean(np.abs(eta) <= .01)))

    if expect_prior == 'uniform' and support > 0:
        d, pvalue = uniform_ks(eta, support)
        stats.update(ks_statistic=d, ks_p_value=pvalue, expected_std=support/math.sqrt(3.))
        check(checks, 'T1 eta is uniform on the ADS domain', pvalue >= alpha,
              f'KS D={d:.5f} p={pvalue:.4f} (alpha={alpha})')
        expected = support/math.sqrt(3.)
        band = 6.*support/math.sqrt(15.*n) if n > 1 else float('inf')
        check(checks, 'T1 eta std matches uniform variance', abs(std-expected) <= band,
              f'std={std:.6f} expected={expected:.6f} band={band:.6f}')
        expected_zero = max(support-.99, 0.)/(2.*support)
        stats.update(expected_fraction_srp_scale_le_1pct=expected_zero)
        zero_band = 5.*math.sqrt(max(expected_zero*(1.-expected_zero), 1e-12)/n)
        check(checks, 'T1 populates the SRP~0 tail at the uniform rate',
              expected_zero > 0. and abs(near_zero-expected_zero) <= zero_band,
              f'P(1+eta<=1%)={near_zero:.5f} expected={expected_zero:.5f} band={zero_band:.5f}')
    elif support > 0:
        stats.update(expected_std=None, expected_fraction_srp_scale_le_1pct=None)
        check(checks, 'T1 eta is not degenerate', np.max(np.abs(eta)) > 0.,
              'legacy normal prior: 1+eta->0 sits at the truncated support edge')

    for name in CASES[:4]:
        if name not in cases:
            continue
        rows = read_points(root/orbit/name/'initial_points.csv')
        eta_c = np.asarray(rows['eta'], dtype=float)
        check(checks, f'{name} keeps eta exactly zero',
              len(eta_c) == 0 or np.all(eta_c == 0.),
              f'max|eta|={np.max(np.abs(eta_c)) if len(eta_c) else 0.:.6e}')

    hashes = {name: cases[name].get('input_state_sha256') for name in cases}
    shared = len(set(hashes.values())) == 1 and hashes['T1'] is not None
    check(checks, 'C0..C3 and T1 share one bit-identical 6D design', shared,
          json.dumps(hashes, sort_keys=True))
    joint = {name: cases[name].get('joint_population_sha256') for name in cases}
    check(checks, 'T1 joint population differs from C0', joint['T1'] != joint['C0'],
          'only the eta column may differ; equal hashes would mean eta was dropped')
    stats['input_state_sha256'] = hashes['T1']
    return dict(orbit=orbit, passed=all(c['ok'] for c in checks),
                checks=checks, statistics=stats)


def parser():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--root', type=Path, required=True, help='campaign output root, e.g. ~/SRP_effect')
    p.add_argument('--orbits', nargs='+', help='orbit subdirectories (default: every orbit with a T1 case)')
    p.add_argument('--expect-prior', choices=('normal', 'uniform'), default='uniform')
    p.add_argument('--eta-sigma', type=float,
                   help='optional cross-check that the recorded support equals domain_sigma*eta_sigma')
    p.add_argument('--domain-sigma', type=float, default=3.)
    p.add_argument('--alpha', type=float, default=.001,
                   help='KS significance level; 1e-3 keeps the false-alarm rate near 1.7%% over 17 orbits')
    p.add_argument('--json', type=Path, help='write the full report here')
    return p


def main(argv=None):
    args = parser().parse_args(argv)
    root = args.root.expanduser().resolve()
    if not root.is_dir():
        raise ValueError(f'campaign root does not exist: {root}')
    orbits = args.orbits or sorted(p.parent.parent.name for p in root.glob('*/T1/population.json'))
    if not orbits:
        raise ValueError(f'no <orbit>/T1/population.json under {root}')
    reports = [check_orbit(root, orbit, args.expect_prior, args.eta_sigma,
                           args.domain_sigma, args.alpha) for orbit in orbits]
    for report in reports:
        stats = report['statistics']
        print(f"[{report['orbit']}] prior={stats.get('eta_prior')} n={stats.get('sample_count')} "
              f"eta=[{stats.get('eta_min', float('nan')):+.6f}, {stats.get('eta_max', float('nan')):+.6f}] "
              f"mean={stats.get('eta_mean', float('nan')):+.2e} "
              f"std={stats.get('eta_std', float('nan')):.6f} "
              f"P(1+eta<=1%)={stats.get('fraction_srp_scale_le_1pct', float('nan')):.4%} "
              f"P(|eta|<=1%)={stats.get('fraction_abs_eta_le_1pct', float('nan')):.4%}")
        for item in report['checks']:
            print(f"  {'PASS' if item['ok'] else 'FAIL'}  {item['name']}: {item['detail']}")
    failed = [r['orbit'] for r in reports if not r['passed']]
    print(f'{len(reports)-len(failed)}/{len(reports)} orbits pass the T1 eta-prior audit')
    if args.json:
        args.json.write_text(json.dumps(dict(schema='t1_eta_prior_audit_v1', root=str(root),
                                             expect_prior=args.expect_prior,
                                             eta_sigma=args.eta_sigma,
                                             domain_sigma=args.domain_sigma, alpha=args.alpha,
                                             orbits=reports), indent=2)+'\n', encoding='utf-8')
    return 0 if not failed else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception as exc:
        print(f'ERROR: {exc}', file=sys.stderr)
        sys.exit(2)
