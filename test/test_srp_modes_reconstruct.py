import csv
from pathlib import Path
import tempfile
import unittest

from test_srp_ads_reconstruct import load_module


def write_csv(path, fields, rows):
    with path.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def fixture(prefix, split=True):
    fields = ['patch_id', 'time_seconds'] + [f'center{i}' for i in range(1,9)] + [f'width{i}' for i in range(1,9)] + ['split_count']
    rows, histories, terms = [], [], []
    coeff_fields = ['patch_id','component'] + [f'power{i}' for i in range(1,9)] + ['coefficient']
    for patch_id, c8 in ([(1,-0.5),(2,0.5)] if split else [(1,0.0)]):
        row = dict(patch_id=patch_id,time_seconds=3600,split_count=int(split))
        row.update({f'center{i}': c8 if i==8 else 0 for i in range(1,9)})
        row.update({f'width{i}': 1 if i==8 and split else 2 for i in range(1,9)})
        rows.append(row)
        if split:
            histories.append(dict(patch_id=patch_id,sequence=1,split_direction=-8 if patch_id==1 else 8))
        def term(component, coefficient, p7=0, p8=0):
            powers = {f'power{i}': 0 for i in range(1,9)}
            powers.update(power7=p7,power8=p8)
            terms.append(dict(patch_id=patch_id,component=component,coefficient=coefficient,**powers))
        span=0.5 if split else 1
        term(1,10+3*c8); term(1,2+c8,1); term(1,3*span,0,1); term(1,span,1,1)
        for component in range(2,7):
            term(component,10*component)
    write_csv(Path(f'{prefix}_meta.csv'),fields,rows)
    write_csv(Path(f'{prefix}_coeff.csv'),coeff_fields,terms)
    write_csv(Path(f'{prefix}_splits.csv'),['patch_id','sequence','split_direction'],histories)
    mapping=[dict(key=k,i=0,j=0,value=v) for k,v in dict(schema_version=2,n_variables=8,da_variables=8,parameterization='SRP_MODES',
        da_order=2,domain_sigma=3,eta_sigma=0,mode_csv_sha256='a'*64,mode_metadata_sha256='b'*64).items()]
    mapping += [dict(key='variable_role',i=i,j=0,value=role) for i,role in enumerate(['x','y','z','vx','vy','vz','s','u'],1)]
    write_csv(Path(f'{prefix}_map.csv'),['key','i','j','value'],mapping)


class ModesSnapshots(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.reconstruct=load_module()

    def test_8d_split_polynomial_and_physical_round_trip(self):
        with tempfile.TemporaryDirectory() as folder:
            prefix=Path(folder)/'modes'
            fixture(prefix)
            snap=self.reconstruct.load_snapshot(prefix)
            self.assertEqual(snap.n_variables,8)
            self.assertEqual(snap.parameterization,'SRP_MODES')
            for q7,q8 in [(-1,-1),(1,1),(0.2,-0.75),(-0.3,0.75),(0,0)]:
                point=(0,)*6+(q7,q8)
                _,state=snap.evaluate_global(point)
                self.assertAlmostEqual(state[0],10+2*q7+3*q8+q7*q8,places=13)
                s,u=snap.physical_modes(q7,q8)
                self.assertAlmostEqual(2*s-1,q7)
                self.assertAlmostEqual(2*u-1,q8)
            samples=self.reconstruct.generate_local_samples(8,512,0.5)
            self.assertEqual(len(samples),512)
            self.assertTrue(all(len(sample.coordinates)==8 for sample in samples))
            result=self.reconstruct.reconstruct_snapshot(prefix,points_per_patch=512,absolute_only=True,verify=False)
            self.assertEqual(result['n_variables'],8)

    def test_bad_v2_contract_is_rejected(self):
        for key,value in [('mode_csv_sha256','bad'),('parameterization','SRP_ETA'),('da_variables','7'),('schema_version','3')]:
            with self.subTest(key=key),tempfile.TemporaryDirectory() as folder:
                prefix=Path(folder)/'modes'; fixture(prefix)
                path=Path(f'{prefix}_map.csv')
                with path.open(newline='') as stream:
                    rows=list(csv.DictReader(stream))
                for row in rows:
                    if row['key']==key: row['value']=value
                write_csv(path,['key','i','j','value'],rows)
                with self.assertRaises(self.reconstruct.SnapshotError):
                    self.reconstruct.load_snapshot(prefix)

    def test_v2_active_six_or_seven_embedded_in_eight(self):
        for active in (6,7):
            with self.subTest(active=active),tempfile.TemporaryDirectory() as folder:
                prefix=Path(folder)/'embedded'; fixture(prefix,split=False)
                path=Path(f'{prefix}_map.csv')
                with path.open(newline='') as stream: rows=list(csv.DictReader(stream))
                for row in rows:
                    if row['key']=='n_variables': row['value']=active
                    if row['key']=='parameterization': row['value']='NONE' if active==6 else 'SRP_ETA'
                    if row['key']=='eta_sigma': row['value']=0 if active==6 else 0.1
                    if row['key']=='variable_role' and int(row['i'])>6:
                        row['value']='eta' if active==7 and int(row['i'])==7 else 'inactive'
                write_csv(path,['key','i','j','value'],rows)
                coeff_path=Path(f'{prefix}_coeff.csv')
                with coeff_path.open(newline='') as stream:
                    reader=csv.DictReader(stream); fields=reader.fieldnames
                    terms=[row for row in reader if not any(int(row[f'power{i}']) for i in range(active+1,9))]
                write_csv(coeff_path,fields,terms)
                output=Path(folder)/'samples.csv'
                self.reconstruct.reconstruct_snapshot(prefix,output=output,points_per_patch=512,absolute_only=True,verify=False)
                with output.open(newline='') as stream:
                    samples=list(csv.DictReader(stream))
                self.assertEqual(len(samples),512)
                self.assertTrue(all(None not in row and all(value is not None for value in row.values()) for row in samples))
                for row in samples:
                    self.assertEqual(float(row['local8']),0)
                    self.assertEqual(float(row['u8']),0)
                    self.assertAlmostEqual(float(row['x_km']),10+2*float(row['u7']) if active==7 else 10)


if __name__=='__main__':
    unittest.main()
