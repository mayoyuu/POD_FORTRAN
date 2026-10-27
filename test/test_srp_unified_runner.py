import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from run_srp_unified import resolve_case, parse_method, read_opm, population_for_case, parser, load_samples
from srp_unified_math import shared_design


class RunnerContractTests(unittest.TestCase):
    def test_build_only_does_not_require_output_directory(self):
        args = parser().parse_args(['--build-only'])
        self.assertTrue(args.build_only)
        self.assertIsNone(args.output)

    def test_case_matrix_is_unambiguous(self):
        self.assertEqual(resolve_case('C0')['model'], 'off')
        self.assertEqual(resolve_case('C1')['model'], 'cannonball')
        self.assertEqual(resolve_case('C2')['attitude'], 'moon')
        self.assertEqual(resolve_case('C3')['attitude'], 'earth')
        self.assertEqual(resolve_case('T1')['error'], 'eta')
        self.assertEqual(resolve_case('T2')['error'], 'attitude')
        with self.assertRaises(ValueError): resolve_case('C8')

    def test_method_selection_preserves_order_and_adaptation(self):
        self.assertEqual(parse_method('DAMC-3', 4), ('DAMC', 3, False))
        self.assertEqual(parse_method('DAMC-auto', 2), ('DAMC', 2, True))
        self.assertEqual(parse_method('UT', 4), ('UT', 4, False))
        with self.assertRaises(ValueError): parse_method('damc--1', 4)

    def test_case_populations_share_initial_state(self):
        design = shared_design(500, 5)
        mean = np.arange(6, dtype=float)
        cov = np.diag([10., 10., 10., 1e-9, 1e-9, 1e-9])
        pops = [population_for_case(design, mean, cov, .3, resolve_case(c))
                for c in ('C0', 'C1', 'C2', 'C3', 'T1', 'T2')]
        for p in pops[1:]: np.testing.assert_array_equal(p[:, :6], pops[0][:, :6])
        np.testing.assert_array_equal(pops[4][:, 6], .3*design[:, 6])
        np.testing.assert_array_equal(pops[5][:, 7:], design[:, 7:])

    def test_uniform_eta_prior_population_is_bounded_and_shares_state(self):
        design, physical, ids = load_samples(None, 4000, 20261021, 3.0, 1/3, 'uniform')
        self.assertIsNone(physical)
        np.testing.assert_array_equal(ids, np.arange(1, 4001))
        pops = [population_for_case(design, np.zeros(6), np.eye(6), 1/3, resolve_case(c))
                for c in ('C0', 'C1', 'C2', 'C3', 'T1')]
        for p in pops[:4]: self.assertTrue(np.all(p[:, 6] == 0))
        np.testing.assert_array_equal(pops[4][:, 6], (1/3)*design[:, 6])
        self.assertLessEqual(np.abs(pops[4][:, 6]).max(), 1+1e-12)
        self.assertGreater(np.abs(pops[4][:, 6]).max(), .99)
        self.assertAlmostEqual(pops[4][:, 6].std(), 1/np.sqrt(3), delta=.02)
        self.assertTrue(.002 < np.mean(1+pops[4][:, 6] <= .01) < .012)
        for p in pops[1:]: np.testing.assert_array_equal(p[:, :6], pops[0][:, :6])

    def test_uniform_eta_prior_outside_physical_support_is_rejected(self):
        with self.assertRaises(ValueError): load_samples(None, 100, 1, 3.0, .5, 'uniform')
        design, _, _ = load_samples(None, 100, 1, 3.0, .3, 'uniform')
        self.assertLessEqual(np.abs(design[:, 6]).max(), 3)

    def test_eta_prior_cli_defaults_to_normal_and_accepts_uniform(self):
        self.assertEqual(parser().parse_args(['--build-only']).eta_prior, 'normal')
        self.assertEqual(parser().parse_args(['--build-only', '--eta-prior', 'uniform']).eta_prior, 'uniform')

    def test_opm_reads_full_covariance_and_validates_spd(self):
        root = Path(__file__).resolve().parents[1]
        mean, cov, data = read_opm(root/'OPM/L1Halo-1/L1Halo-1_init.opm.json', 2.)
        self.assertEqual(mean.shape, (6,))
        self.assertAlmostEqual(np.trace(cov[:3, :3]), 400.)
        self.assertEqual(data['ASL_CAT_ID'], 'L1Halo-1')



class CompilerBuildDiscoveryTests(unittest.TestCase):
    def check_selected_build(self, selected_mtime, unrelated_mtime):
        import os
        from types import SimpleNamespace
        from unittest.mock import patch
        import run_srp_unified as runner

        with tempfile.TemporaryDirectory(prefix='srp build ') as tmp:
            root = Path(tmp)
            unrelated = root/'build/gfortran_old/app/run_srp_unified'
            selected = root/'build/gfortran-13_current/app/run_srp_unified'
            for binary in (unrelated, selected):
                binary.parent.mkdir(parents=True)
                binary.write_text('test executable')
                binary.chmod(0o755)
            os.utime(unrelated, (unrelated_mtime, unrelated_mtime))
            os.utime(selected, (selected_mtime, selected_mtime))

            def fake_fpm(*args, **kwargs):
                kwargs['stdout'].write('__POD_UNIFIED_TARGET__ '+str(selected)+'\n')
                return SimpleNamespace(returncode=0)

            with patch.object(runner, 'ROOT', root), patch.object(runner.subprocess, 'run', fake_fpm):
                executable, _ = runner.build_executable(SimpleNamespace(executable=None), {})
            self.assertEqual(executable, selected.resolve())

    def test_build_finds_named_gfortran_compiler_instead_of_stale_binary(self):
        self.check_selected_build(selected_mtime=200, unrelated_mtime=100)

    def test_cached_selected_build_wins_over_newer_unrelated_binary(self):
        self.check_selected_build(selected_mtime=100, unrelated_mtime=200)


if __name__ == '__main__': unittest.main()
