import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from run_srp_unified import resolve_case, parse_method, read_opm, population_for_case, parser
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

    def test_opm_reads_full_covariance_and_validates_spd(self):
        root = Path(__file__).resolve().parents[1]
        mean, cov, data = read_opm(root/'OPM/L1Halo-1/L1Halo-1_init.opm.json', 2.)
        self.assertEqual(mean.shape, (6,))
        self.assertAlmostEqual(np.trace(cov[:3, :3]), 400.)
        self.assertEqual(data['ASL_CAT_ID'], 'L1Halo-1')


if __name__ == '__main__': unittest.main()
