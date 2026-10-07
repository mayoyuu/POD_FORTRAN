#!/usr/bin/env python3
"""Pack SRP paired point staging CSV into a validated, compressed HDF5 file."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import tempfile
import h5py
import numpy as np

SOURCE_SCHEMA='srp_mc_points_v1'
OUTPUT_SCHEMA='srp_mc_points_hdf5_v1'
COMPONENTS=('dx_km','dy_km','dz_km','dvx_kms','dvy_kms','dvz_kms')
VALUE_COLUMNS=tuple(f'{group}_{component}'
                    for group in ('treatment','control','delta')
                    for component in COMPONENTS)


class PointFormatError(ValueError):
    pass


def _sha256(path: Path) -> str:
    digest=hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def _rows(path: Path):
    handle=path.open(newline='')
    try:
        schema=handle.readline().rstrip('\r\n')
        if schema!=f'schema={SOURCE_SCHEMA}':
            raise PointFormatError('point CSV schema line is missing or unsupported')
        reader=csv.DictReader(handle)
        required={'hour','keyframe_kind','sample_id',*VALUE_COLUMNS}
        if reader.fieldnames is None or set(reader.fieldnames)!=required:
            raise PointFormatError('point CSV header does not match srp_mc_points_v1')
        yield from reader
    finally:
        handle.close()


def _parse_row(row: dict[str,str],row_number: int):
    try:
        hour=float(row['hour'])
        sample_id=int(row['sample_id'])
        kind=row['keyframe_kind']
        values=np.asarray([float(row[name]) for name in VALUE_COLUMNS],dtype=np.float64)
    except (KeyError,ValueError,TypeError) as exc:
        raise PointFormatError(f'malformed point row {row_number}') from exc
    if not math.isfinite(hour) or sample_id<1 or not kind or ',' in kind:
        raise PointFormatError(f'invalid point identity at row {row_number}')
    if not np.isfinite(values).all():
        raise PointFormatError(f'nonfinite point value at row {row_number}')
    treatment=values[:6]
    control=values[6:12]
    delta=values[12:]
    tolerance=1e-12*np.maximum(1.0,np.maximum(np.abs(delta),
                                             np.abs(treatment-control)))
    if np.any(np.abs(delta-(treatment-control))>tolerance):
        raise PointFormatError(f'paired delta mismatch at row {row_number}')
    return hour,kind,sample_id,treatment,control,delta


def _scan(path: Path,expected_per_frame: int):
    if expected_per_frame<1:
        raise PointFormatError('expected point count must be positive')
    frames=[]
    frame_keys=set()
    current_key=None
    current_ids=[]
    baseline_ids=None
    row_count=0

    def finish_frame():
        nonlocal baseline_ids,current_ids,current_key
        if current_key is None:
            return
        if len(current_ids)!=expected_per_frame:
            raise PointFormatError(
                f'frame {current_key!r} has {len(current_ids)} points; '
                f'expected {expected_per_frame}')
        if len(set(current_ids))!=len(current_ids):
            raise PointFormatError(f'duplicate sample_id in frame {current_key!r}')
        if baseline_ids is None:
            baseline_ids=tuple(current_ids)
        elif tuple(current_ids)!=baseline_ids:
            raise PointFormatError(f'sample_id order changed in frame {current_key!r}')
        frames.append(current_key)

    for row_count,row in enumerate(_rows(path),start=1):
        hour,kind,sample_id,*_= _parse_row(row,row_count)
        key=(hour,kind)
        if current_key is None:
            current_key=key
        elif key!=current_key:
            finish_frame()
            frame_keys.add(current_key)
            if key in frame_keys:
                raise PointFormatError(f'noncontiguous duplicate frame {key!r}')
            current_key=key
            current_ids=[]
        current_ids.append(sample_id)
    finish_frame()
    if row_count==0:
        raise PointFormatError('point CSV contains no data rows')
    if len(set(frames))!=len(frames):
        raise PointFormatError('duplicate frame key')
    return frames,np.asarray(baseline_ids,dtype=np.int64),row_count


def pack_points(source: str|Path,output: str|Path,expected_per_frame: int,
                point_role: str,report_json: str|Path|None=None) -> dict:
    source=Path(source)
    output=Path(output)
    if point_role not in ('visualization','keyframe'):
        raise PointFormatError('point role must be visualization or keyframe')
    frames,sample_ids,row_count=_scan(source,expected_per_frame)
    output.parent.mkdir(parents=True,exist_ok=True)
    fd,temp_name=tempfile.mkstemp(prefix=output.name+'.',suffix='.tmp',
                                  dir=output.parent)
    os.close(fd)
    temp_path=Path(temp_name)
    frame_count=len(frames)
    chunks=(1,min(expected_per_frame,2048),6)
    try:
        with h5py.File(temp_path,'w') as h5:
            h5.attrs['output_schema_version']=OUTPUT_SCHEMA
            h5.attrs['source_schema_version']=SOURCE_SCHEMA
            h5.attrs['point_role']=point_role
            h5.attrs['source_csv_sha256']=_sha256(source)
            h5.attrs['source_csv_name']=source.name
            h5.attrs['frame_count']=frame_count
            h5.attrs['point_count_per_frame']=expected_per_frame
            h5.attrs['state_units']='km,km,km,km/s,km/s,km/s'
            if report_json is not None:
                report_path=Path(report_json)
                report=json.loads(report_path.read_text())
                h5.attrs['paired_report_schema']=str(
                    report.get('output_schema_version',''))
                h5.attrs['paired_report_sha256']=_sha256(report_path)
            h5.create_dataset('hour',data=np.asarray([key[0] for key in frames]),
                              dtype='f8')
            string_dtype=h5py.string_dtype(encoding='utf-8')
            h5.create_dataset('keyframe_kind',
                              data=np.asarray([key[1] for key in frames],
                                              dtype=object),
                              dtype=string_dtype)
            id_dataset=h5.create_dataset(
                'sample_id',shape=(frame_count,expected_per_frame),dtype='i8',
                chunks=(1,min(expected_per_frame,8192)),compression='gzip',
                compression_opts=4,shuffle=True)
            treatment_dataset=h5.create_dataset(
                'treatment_error',shape=(frame_count,expected_per_frame,6),
                dtype='f8',chunks=chunks,compression='gzip',
                compression_opts=4,shuffle=True)
            control_dataset=h5.create_dataset(
                'control_error',shape=(frame_count,expected_per_frame,6),
                dtype='f8',chunks=chunks,compression='gzip',
                compression_opts=4,shuffle=True)
            delta_dataset=h5.create_dataset(
                'delta_error',shape=(frame_count,expected_per_frame,6),
                dtype='f8',chunks=chunks,compression='gzip',
                compression_opts=4,shuffle=True)

            treatment=np.empty((expected_per_frame,6),dtype=np.float64)
            control=np.empty_like(treatment)
            delta=np.empty_like(treatment)
            frame_index=0
            point_index=0
            active_key=frames[0]
            for row_number,row in enumerate(_rows(source),start=1):
                hour,kind,sample_id,t,c,d=_parse_row(row,row_number)
                key=(hour,kind)
                if key!=active_key:
                    id_dataset[frame_index,:]=sample_ids
                    treatment_dataset[frame_index,:,:]=treatment
                    control_dataset[frame_index,:,:]=control
                    delta_dataset[frame_index,:,:]=delta
                    frame_index+=1
                    active_key=key
                    point_index=0
                if sample_id!=sample_ids[point_index]:
                    raise PointFormatError('sample order changed between validation passes')
                treatment[point_index,:]=t
                control[point_index,:]=c
                delta[point_index,:]=d
                point_index+=1
            id_dataset[frame_index,:]=sample_ids
            treatment_dataset[frame_index,:,:]=treatment
            control_dataset[frame_index,:,:]=control
            delta_dataset[frame_index,:,:]=delta
            h5.flush()
        os.replace(temp_path,output)
    except BaseException:
        try: temp_path.unlink()
        except FileNotFoundError: pass
        raise
    return {'output_schema_version':OUTPUT_SCHEMA,'source':str(source),
            'output':str(output),'source_csv_sha256':_sha256(source),
            'frame_count':frame_count,'point_count':expected_per_frame,
            'row_count':row_count}


def main() -> int:
    parser=argparse.ArgumentParser()
    parser.add_argument('--input',required=True,type=Path)
    parser.add_argument('--output',required=True,type=Path)
    parser.add_argument('--expected-per-frame',required=True,type=int)
    parser.add_argument('--point-role',required=True,
                        choices=('visualization','keyframe'))
    parser.add_argument('--report-json',type=Path)
    parser.add_argument('--manifest',type=Path)
    args=parser.parse_args()
    result=pack_points(args.input,args.output,args.expected_per_frame,
                       args.point_role,args.report_json)
    text=json.dumps(result,indent=2)+'\n'
    if args.manifest:
        args.manifest.write_text(text)
    print(text,end='')
    return 0


if __name__=='__main__':
    raise SystemExit(main())
