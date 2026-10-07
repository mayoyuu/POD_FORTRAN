#!/usr/bin/env python3
import csv
import tempfile
import unittest
from pathlib import Path
import sys
import h5py
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
import pack_srp_mc_points_hdf5 as packer

HEADER=['hour','keyframe_kind','sample_id',
        'treatment_dx_km','treatment_dy_km','treatment_dz_km',
        'treatment_dvx_kms','treatment_dvy_kms','treatment_dvz_kms',
        'control_dx_km','control_dy_km','control_dz_km',
        'control_dvx_kms','control_dvy_kms','control_dvz_kms',
        'delta_dx_km','delta_dy_km','delta_dz_km',
        'delta_dvx_kms','delta_dvy_kms','delta_dvz_kms']


class Hdf5PackerTest(unittest.TestCase):
    def make_csv(self,path,duplicate=False,bad_delta=False):
        with Path(path).open('w',newline='') as f:
            f.write('schema=srp_mc_points_v1\n')
            writer=csv.DictWriter(f,fieldnames=HEADER)
            writer.writeheader()
            for hour in (0,1):
                for sample in (2,5,9):
                    treatment=np.arange(6,dtype=float)+hour+sample
                    control=treatment-0.25
                    delta=treatment-control
                    if bad_delta and hour==1 and sample==5:
                        delta[0]=99
                    row={'hour':hour,'keyframe_kind':'hourly_visualization',
                         'sample_id':sample}
                    for prefix,values in (('treatment',treatment),('control',control),
                                          ('delta',delta)):
                        for name,value in zip(('dx_km','dy_km','dz_km',
                                               'dvx_kms','dvy_kms','dvz_kms'),values):
                            row[f'{prefix}_{name}']=value
                    writer.writerow(row)
                    if duplicate and hour==0 and sample==2:
                        writer.writerow(row)

    def test_pack_preserves_pairing_and_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'points.csv'
            output=Path(tmp)/'points.h5'
            self.make_csv(source)
            result=packer.pack_points(source,output,3,'visualization')
            self.assertEqual(result['frame_count'],2)
            self.assertEqual(result['point_count'],3)
            with h5py.File(output,'r') as h5:
                self.assertEqual(h5.attrs['output_schema_version'],
                                 'srp_mc_points_hdf5_v1')
                self.assertEqual(h5.attrs['point_role'],'visualization')
                np.testing.assert_array_equal(h5['sample_id'][0],[2,5,9])
                np.testing.assert_allclose(h5['delta_error'][:],
                                           h5['treatment_error'][:]-h5['control_error'][:])
                self.assertEqual(h5['treatment_error'].shape,(2,3,6))

    def test_duplicate_point_is_rejected_atomically(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'points.csv'
            output=Path(tmp)/'points.h5'
            self.make_csv(source,duplicate=True)
            with self.assertRaises(packer.PointFormatError):
                packer.pack_points(source,output,3,'visualization')
            self.assertFalse(output.exists())

    def test_inconsistent_delta_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            source=Path(tmp)/'points.csv'
            output=Path(tmp)/'points.h5'
            self.make_csv(source,bad_delta=True)
            with self.assertRaises(packer.PointFormatError):
                packer.pack_points(source,output,3,'visualization')
            self.assertFalse(output.exists())


if __name__=='__main__':
    unittest.main()
