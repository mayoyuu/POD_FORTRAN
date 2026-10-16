#!/usr/bin/env python3
"""Read-only P4 case audit: numerical validity, convergence, and proxy-model error."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
STATE=('x_km','y_km','z_km','vx_kms','vy_kms','vz_kms')
NOMINAL=tuple('nominal_'+k for k in STATE)
ERROR=('dx_km','dy_km','dz_km','dvx_kms','dvy_kms','dvz_kms')


class CaseError(RuntimeError): pass


def require(condition,message):
    if not condition: raise CaseError(message)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def load_json(path):
    return json.loads(Path(path).read_text(),parse_constant=lambda value: (_ for _ in ()).throw(CaseError('nonfinite JSON')))


def rows(path):
    with Path(path).open(newline='') as stream:
        reader=csv.DictReader(stream)
        for row in reader:
            require(None not in row and all(v is not None for v in row.values()),f'invalid CSV fields: {path}')
            yield row


def finite(value):
    result=float(value); require(math.isfinite(result),'nonfinite number'); return result


def vector(row,keys): return tuple(finite(row[k]) for k in keys)


def audit_manifest(folder,report,required):
    manifest=report.get('artifacts_sha256',{})
    require(set(required)<=manifest.keys(),'required artifacts missing from hash manifest')
    for name,expected in manifest.items():
        require(Path(name).name==name,'artifact path must stay in result directory')
        require(sha(folder/name)==expected,f'artifact hash mismatch: {name}')
    for key,item in report['provenance']['inputs'].items():
        require(sha(item['path'])==item['sha256'],f'input hash mismatch: {key}')


def load_run(folder):
    folder=Path(folder).resolve(); report=load_json(folder/'c2_report.json')
    require(report.get('schema_version')=='srp_modes_ads_v1' and report.get('formal_result') is True,'not a formal C2 result')
    require(report.get('n_variables')==8 and report.get('parameterization')=='SRP_MODES' and
            report.get('variable_roles')==['x','y','z','vx','vy','vz','s','u'],'C2 dimension/roles mismatch')
    require(report.get('validation_failed') is False and report.get('precision_unconfirmed') is False,'failed numerical status')
    require(report.get('da_handles_baseline')==report.get('da_handles_final')==0,'DA handle leak')
    require(report.get('depth_limited_patches')==0,'depth-limited result')
    for key,ptol,vtol in (('ads_numerical_error',.1,1e-6),('zero_mode_error',1e-6,1e-9)):
        require(0<=finite(report[key]['max_position_km'])<=ptol and
                0<=finite(report[key]['max_velocity_kms'])<=vtol,f'{key} gate failed')
    audit_manifest(folder,report,('c2_history.csv','c2_shape.csv','c2_initial_points.csv','c2_probes.csv'))
    cadence=finite(report['save_hours']); duration=finite(report['duration_hours'])
    require(cadence>0 and duration>0 and math.isclose(round(duration/cadence)*cadence,duration,abs_tol=1e-10),'invalid hour grid')
    hours={i*cadence for i in range(round(duration/cadence)+1)}
    nominal={}
    for row in rows(folder/'c2_history.csv'):
        hour=finite(row['hour']); require(hour in hours and hour not in nominal,'missing/duplicate/off-grid history hour')
        require(math.isclose(finite(row['epoch_et']),report['epoch0_tdb_s']+hour*3600,rel_tol=0,abs_tol=1e-6),'history epoch mismatch')
        nominal[hour]=vector(row,NOMINAL)
    require(nominal.keys()==hours,'incomplete history hours')
    points={}
    for row in rows(folder/'c2_initial_points.csv'):
        index=int(row['point_id']); require(index not in points and index>0,'duplicate/invalid initial point')
        points[index]=vector(row,[f'u{i}' for i in range(1,9)])
        require(all(abs(q)<=1 for q in points[index]),'point outside unit domain')
    require(len(points)==report['shape_count'],'incomplete initial point grid')
    shape={hour:{} for hour in hours}
    for row in rows(folder/'c2_shape.csv'):
        hour=finite(row['hour']); point=int(row['point_id'])
        require(hour in hours and point in points and point not in shape[hour],'duplicate/off-grid shape point')
        shape[hour][point]=vector(row,ERROR)
    require(all(cloud.keys()==points.keys() for cloud in shape.values()),'incomplete shape hours/points')
    probes={hour:set() for hour in hours}; maxima=[0.,0.]
    for row in rows(folder/'c2_probes.csv'):
        hour=finite(row['hour']); point=int(row['point_id'])
        require(hour in hours and point in points and point not in probes[hour],'duplicate/off-grid probe')
        probes[hour].add(point)
        error=norms(difference(vector(row,['ads_'+k for k in ('x','y','z','vx','vy','vz')]),
                               vector(row,['real_'+k for k in ('x','y','z','vx','vy','vz')])))
        require(error[0]<=.1 and error[1]<=1e-6,'actual probe numerical gate failed')
        maxima=[max(a,b) for a,b in zip(maxima,error)]
    require(all(len(ids)==report['validation_count'] and ids==probes[min(hours)] for ids in probes.values()),
            'incomplete/inconsistent probe grid')
    require(math.isclose(maxima[0],report['ads_numerical_error']['max_position_km'],rel_tol=1e-6,abs_tol=1e-10) and
            math.isclose(maxima[1],report['ads_numerical_error']['max_velocity_kms'],rel_tol=1e-6,abs_tol=1e-14),
            'probe maxima/report disagreement')
    return dict(folder=folder,report=report,nominal=nominal,shape=shape,points=points,probe_ids=probes[min(hours)])


def difference(a,b): return tuple(x-y for x,y in zip(a,b))
def norms(state): return math.hypot(*state[:3]),math.hypot(*state[3:])


def envelopes(run,reference):
    return {hour:tuple(max(norms(tuple(offset[i]+run['nominal'][hour][i]-reference[hour][i]
                         for i in range(6)))[component] for offset in cloud.values())
                         for component in (0,1)) for hour,cloud in run['shape'].items()}


def vector_metrics(left,right):
    require(left.keys()==right.keys(),'trajectory hour mismatch')
    values={hour:norms(difference(left[hour],right[hour])) for hour in left}
    ph=max(values,key=lambda hour:values[hour][0]); vh=max(values,key=lambda hour:values[hour][1])
    final=max(values)
    return dict(max_position_km=values[ph][0],max_position_hour=ph,max_velocity_kms=values[vh][1],
                max_velocity_hour=vh,final_position_km=values[final][0],final_velocity_kms=values[final][1])


def compare_runs(baseline,other,kind,reference=None):
    require(kind in ('same','half_hour','tight_rkf','da5'),'unknown comparison type')
    a,b=baseline['report'],other['report']
    for key in ('n_variables','parameterization','variable_roles','epoch0_tdb_s','duration_hours','save_hours',
                'shape_count','validation_count','domain_sigma','max_depth','position_tolerance_km',
                'velocity_tolerance_kms','dt_min_s','dt_max_s'):
        require(a[key]==b[key],f'run identity mismatch: {key}')
    require(baseline['points']==other['points'],'study point grid mismatch')
    require(baseline['probe_ids']==other['probe_ids'],'Real validation point grid mismatch')
    ai,bi=a['provenance']['inputs'],b['provenance']['inputs']
    for key in ('opm','mode_csv','mode_metadata'):
        if key!='opm' and kind=='half_hour': continue
        require(ai[key]['sha256']==bi[key]['sha256'],f'input identity mismatch: {key}')
    ac,bc=a['provenance']['source_files_sha256'],b['provenance']['source_files_sha256']
    require({k:v for k,v in ac.items() if k.startswith('src/')}==
            {k:v for k,v in bc.items() if k.startswith('src/')},'core source identity mismatch')
    require(baseline['nominal'].keys()==other['nominal'].keys(),'comparison hours mismatch')
    reference=baseline['nominal'] if reference is None else reference
    ea=envelopes(baseline,reference); eb=envelopes(other,reference)
    changes=[]
    for component in (0,1):
        scale=max(value[component] for value in ea.values())
        delta=max(abs(ea[h][component]-eb[h][component]) for h in ea)
        changes.append(delta/scale if scale else (0 if delta==0 else math.inf))
    point_differences=[norms(tuple(other['shape'][hour][p][i]+other['nominal'][hour][i]-
                       baseline['shape'][hour][p][i]-baseline['nominal'][hour][i] for i in range(6)))
                       for hour in ea for p in baseline['points']]
    offsets={key:abs(b['conditional_srp_offset'][key]-value)/abs(value) if value else
                    (0 if b['conditional_srp_offset'][key]==0 else math.inf)
             for key,value in a['conditional_srp_offset'].items()}
    return dict(position_envelope_relative_change=changes[0],velocity_envelope_relative_change=changes[1],
        conditional_srp_offset_relative_changes=offsets,
        max_state_difference_position_km=max(x[0] for x in point_differences),
        max_state_difference_velocity_kms=max(x[1] for x in point_differences),
        convergence_confirmed=all(x<=.05 for x in [*changes,*offsets.values()]),threshold=.05)


def import_tool(name):
    spec=importlib.util.spec_from_file_location('case_'+name,ROOT/'tools'/f'{name}.py')
    result=importlib.util.module_from_spec(spec); spec.loader.exec_module(result); return result


def audit_snapshot_points(run):
    """Re-evaluate immutable snapshot polynomials in memory, never rewrite their audit CSVs."""
    tool=import_tool('reconstruct_srp_ads_snapshot'); report=run['report']; audits=[]; observed=set()
    for prefix in report['snapshot_prefixes']:
        snapshot=tool.load_snapshot(run['folder']/prefix); hour=snapshot.time_seconds/3600
        require(hour in run['nominal'] and hour not in observed,'invalid/duplicate snapshot hour')
        require(snapshot.n_variables==8 and snapshot.da_order==report['da_order'] and
                snapshot.domain_sigma==report['domain_sigma'],'snapshot run settings mismatch')
        require(snapshot.mode_hashes=={key:report[key] for key in ('mode_csv_sha256','mode_metadata_sha256')},
                'snapshot mode identity mismatch')
        maxima=[0.,0.]
        for point,coordinates in run['points'].items():
            _,state=snapshot.evaluate_global(coordinates)
            stored=tuple(a+b for a,b in zip(run['shape'][hour][point],run['nominal'][hour]))
            error=norms(difference(state,stored)); maxima=[max(a,b) for a,b in zip(maxima,error)]
        require(maxima[0]<=1e-8 and maxima[1]<=1e-11,'snapshot/saved cloud mismatch')
        audits.append(dict(hour=hour,point_count=len(run['points']),max_position_km=maxima[0],max_velocity_kms=maxima[1]))
        observed.add(hour)
    require({0,report['duration_hours']}<=observed,'initial/final snapshot missing')
    return audits


def audit_controls(c0,c1,baseline):
    """Reuse historical domains; never call a propagator or rewrite their artifacts."""
    c0=Path(c0).resolve(); c1=Path(c1).resolve()
    common=import_tool('run_srp_modes_ads'); snapshots=import_tool('reconstruct_srp_ads_snapshot')
    base_report=baseline['report']; prefix='DRO-2_init'
    c0r=load_json(c0/f'{prefix}_report.json'); c1r=load_json(c1/f'{prefix}_report.json')
    require(c0r['opm_sha256']==base_report['provenance']['inputs']['opm']['sha256'],'C0 OPM identity mismatch')
    require(sha(c0r['source_config'])==c0r['config_sha256'],'C0 preserved config identity mismatch')
    bc=common.read_config(base_report['provenance']['inputs']['config']['path'])
    cc=common.read_config(c0r['source_config'])
    require({k:v for k,v in bc.items() if not k.startswith('srp_')}==
            {k:v for k,v in cc.items() if not k.startswith('srp_')},'C0 non-SRP physics mismatch')
    pressure_identity=None
    for key in ('srp_cannonball_cr','srp_cannonball_area_m2','srp_mass_kg','srp_pressure_1au_n_m2'):
        left,right=finite(bc[key]),finite(cc[key])
        # Config -1 is the existing force model's 1367/c default, not negative pressure.
        if key=='srp_pressure_1au_n_m2':
            if left==-1: left=1367/299792458
            if right==-1: right=1367/299792458
            require(left>0 and right>0,'invalid C0/C2 physical pressure')
            pressure_identity=dict(C2_pressure_n_m2=left,C0_pressure_n_m2=right,
                relative_difference=(left-right)/right,exactly_equal=left==right,
                interpretation='Recorded physical configuration difference; not silently rounded into identity. C0 is a shared comparison reference, not a newly matched control rerun.')
        else:
            require(left==right,f'C0 cannonball physics mismatch: {key}')
    mapping=lambda path:{(row['key'],row['i'],row['j']):finite(row['value']) for row in rows(path)
                         if row['key'] in ('mean','basis','domain_sigma')}
    expected=mapping(baseline['folder']/f'{base_report["snapshot_prefixes"][0]}_map.csv')
    nominal={}; result={}; all_hashes={}
    for name,folder,report in (('C0',c0,c0r),('C1_H',c1,c1r)):
        require(report['duration_hours']==360 and report['save_hours']==1 and report['first_failed_hour'] is None and
                report['first_unconfirmed_hour'] is None,'historical control incomplete/failed')
        require(mapping(folder/f'{prefix}_patch_s0000000_map.csv')==expected,'historical initial mean/basis mismatch')
        initial=snapshots.load_snapshot(folder/f'{prefix}_patch_s0000000')
        final=snapshots.load_snapshot(folder/f'{prefix}_patch_s1296000')
        require(initial.time_seconds==0 and final.time_seconds==1296000,'historical snapshot times mismatch')
        require(final.n_variables==(6 if name=='C0' else 7),'historical dimension mismatch')
        history={}
        for row in rows(folder/f'{prefix}_history.csv'):
            hour=finite(row['hour']); require(hour in baseline['nominal'] and hour not in history,'historical hour mismatch')
            require(math.isclose(finite(row['epoch_et']),base_report['epoch0_tdb_s']+hour*3600,rel_tol=0,abs_tol=1e-6),
                    'historical epoch mismatch')
            require(int(row['validation_failed'])==int(row['precision_unconfirmed'])==int(row['depth_limited'])==0 and
                    finite(row['max_validation_pos_km'])<=.1 and finite(row['max_validation_vel_kms'])<=1e-6,
                    'historical hourly numerical gate failed')
            history[hour]=vector(row,NOMINAL)
        require(history.keys()==baseline['nominal'].keys(),'historical missing hours')
        if name=='C0': nominal=history
        all_hashes[name]={str(path):sha(path) for path in sorted(folder.glob(prefix+'*')) if path.is_file()}
        result[name]=dict(duration_hours=360,hour_count=len(history),source_directory=str(folder),
            historical_numerical_status=report['accuracy_status'],snapshot_n_variables=final.n_variables)
    result['C0']['identity_status']='matched_opm_initial_basis_non_srp_config'
    result['C0']['pressure_comparison']=pressure_identity
    result['C0']['nominal_difference_from_C2']=vector_metrics(nominal,baseline['nominal'])
    finals={name:snapshots.load_snapshot(folder/f'{prefix}_patch_s1296000') for name,folder in (('C0',c0),('C1_H',c1))}
    span=finals['C1_H'].domain_sigma*finals['C1_H'].eta_sigma
    require(math.isclose(span,.9,abs_tol=1e-12),'expected historical C1-H eta +/-0.9')
    envelope={}
    for name in ('C0','C1_30_subdomain','C1_H'):
        source=finals['C0'] if name=='C0' else finals['C1_H']; values=[]
        for point in baseline['points'].values():
            coordinates=point[:6] if name=='C0' else point[:6]+(point[6]*(.3/span if name=='C1_30_subdomain' else 1),)
            _,state=source.evaluate_global(coordinates)
            values.append(norms(difference(state,nominal[360])))
        envelope[name]={'position_km':max(v[0] for v in values),'velocity_kms':max(v[1] for v in values)}
    _,eta_state=finals['C1_H'].evaluate_global((0,)*6+(.3/span,))
    result['C1_30_subdomain']=dict(source_eta_half_span=span,restricted_eta_half_span=.3,
        fully_provenanced_matched_campaign=False,reconstruction_hour=360,
        zero_state_eta_plus_03_position_km=norms(difference(eta_state,nominal[360]))[0],
        zero_state_eta_plus_03_state=list(eta_state),
        limitation='Historical C1-H lacks OPM/config/executable provenance hashes. Matching mean/basis/epoch supports conditional final-time restriction only; no new C1-30 or hourly restricted polynomials are claimed.')
    result['final_matched_state_grid_envelopes']=envelope
    result['input_files_sha256']=all_hashes
    return result,nominal


def analyze(baseline,half_hour,tight_rkf,da5,truth,c0,c1,output):
    output=Path(output).absolute(); require(not output.exists(),'analysis output exists')
    paths=dict(baseline=baseline,half_hour=half_hour,tight_rkf=tight_rkf,da5=da5)
    runs={key:load_run(path) for key,path in paths.items()}; base=runs['baseline']; r=base['report']
    snapshot_audits={key:audit_snapshot_points(run) for key,run in runs.items()}
    require(r['duration_hours']==360 and r['save_hours']==1 and r['shape_count']==512 and r['validation_count']==337,
            'not a complete 360-hour, 512/337-point case')
    require(r['da_order']==4 and r['rel_tol']==r['abs_tol']==1e-12,'unexpected baseline settings')
    common=import_tool('run_srp_modes_ads')
    controls,reference=audit_controls(c0,c1,base)
    cfg=common.read_config(r['provenance']['inputs']['config']['path'])
    comparisons={}
    for kind in ('half_hour','tight_rkf','da5'):
        run=runs[kind]; report=run['report']
        comparisons[kind]=compare_runs(base,run,kind,reference)
        other=common.read_config(report['provenance']['inputs']['config']['path'])
        omit={'srp_mode_table_file','srp_mode_metadata_file'}
        require({k:v for k,v in cfg.items() if k not in omit}=={k:v for k,v in other.items() if k not in omit},
                f'non-mode configuration mismatch: {kind}')
        require(report['da_order']==(5 if kind=='da5' else 4),'DA variant setting mismatch')
        require(report['rel_tol']==report['abs_tol']==(1e-13 if kind=='tight_rkf' else 1e-12),'RKF variant setting mismatch')
        mode_audit,_,_=import_tool('run_srp_modes_truth').preflight(
            report['provenance']['inputs']['opm']['path'],report['provenance']['inputs']['config']['path'])
        require(mode_audit['sample_hours']==(.5 if kind=='half_hour' else 1),'mode-table cadence mismatch')
    truth=Path(truth).resolve(); tr=load_json(truth/'truth_report.json')
    require(tr.get('formal_result') is True,'truth is not formal')
    audit_manifest(truth,tr,('truth_states.csv','runner'))
    require(tr['provenance']['inputs']==r['provenance']['inputs'],'truth/base input identity mismatch')
    truth_tool=import_tool('run_srp_modes_truth')
    initial=[common.load_json(r['provenance']['inputs']['opm']['path'])[k] for k in ('X','Y','Z','X_DOT','Y_DOT','Z_DOT')]
    truth_data=truth_tool.validate_truth_output(truth,tr,360,1,r['epoch0_tdb_s'],initial)
    tracks={scenario:{} for scenario in truth_tool.SCENARIOS}
    for row in truth_data: tracks[row['scenario']][float(row['hour'])]=vector(row,STATE)
    proxy={attitude:vector_metrics(tracks['proxy_'+attitude],tracks['boxwing_'+attitude]) for attitude in ('earth','moon')}
    offsets={name:vector_metrics(track,reference) for name,track in tracks.items()}
    controls['C1_30_subdomain']['independent_Deta_final_error']=dict(zip(('position_km','velocity_kms'),
        norms(difference(controls['C1_30_subdomain']['zero_state_eta_plus_03_state'],tracks['Deta'][360]))))
    snapshot_tool=import_tool('reconstruct_srp_ads_snapshot')
    snapshot=snapshot_tool.load_snapshot(base['folder']/r['snapshot_prefixes'][-1])
    require(snapshot.time_seconds==360*3600,'final snapshot incomplete')
    endpoint={}
    for name,q7,q8 in (('proxy_earth',1,0),('proxy_moon',-1,1)):
        _,state=snapshot.evaluate_global((0,)*6+(q7,q8))
        pos,vel=norms(difference(state,tracks[name][360]))
        require(pos<=.1 and vel<=1e-6,'ADS/independent proxy endpoint numerical gate failed')
        endpoint[name]=dict(position_km=pos,velocity_kms=vel)
    per_hour=[]; eb=envelopes(base,reference)
    for hour in sorted(base['nominal']):
        row={'hour':hour,'baseline_envelope_position_km':eb[hour][0],'baseline_envelope_velocity_kms':eb[hour][1]}
        for attitude in ('earth','moon'):
            pos,vel=norms(difference(tracks['proxy_'+attitude][hour],tracks['boxwing_'+attitude][hour]))
            row[f'{attitude}_proxy_error_position_km']=pos; row[f'{attitude}_proxy_error_velocity_kms']=vel
        per_hour.append(row)
    input_reports={key:{'path':str(run['folder']/'c2_report.json'),'sha256':sha(run['folder']/'c2_report.json')}
                   for key,run in runs.items()}
    input_reports['truth']={'path':str(truth/'truth_report.json'),'sha256':sha(truth/'truth_report.json')}
    result=dict(schema_version='srp_modes_case_analysis_v1',complete_orbit_case=True,duration_hours=360,hour_count=361,
        numerical_gates_passed=True,convergence_confirmed=all(x['convergence_confirmed'] for x in comparisons.values()),
        convergence=comparisons,ads_numerical_error={key:run['report']['ads_numerical_error'] for key,run in runs.items()},
        zero_mode_error={key:run['report']['zero_mode_error'] for key,run in runs.items()},
        surrogate_model_error={'classification':'proxy model limitation, not an ADS numerical gate','attitudes':proxy},
        conditional_srp_offset=offsets,independent_proxy_endpoint_ads_error=endpoint,
        geometry_fallback_count=tr['geometry_fallback_count'],input_reports=input_reports,
        independent_snapshot_audits=snapshot_audits,
        reused_controls=controls,common_reference='reused C0 cannonball nominal hourly history',
        analysis_source_sha256=sha(__file__),interpretation='Sampled deterministic boxes, not a probability distribution or exhaustive continuous-domain guarantee.')
    # Output is exclusively created only after every input and numerical check succeeds.
    output.mkdir(parents=True,exist_ok=False)
    with (output/'hourly_metrics.csv').open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(per_hour[0])); writer.writeheader(); writer.writerows(per_hour)
    result['artifacts_sha256']={'hourly_metrics.csv':sha(output/'hourly_metrics.csv')}
    (output/'case_report.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('baseline','half-hour','tight-rkf','da5','truth','c0','c1','output'): parser.add_argument('--'+key,type=Path,required=True)
    try: report=analyze(**vars(parser.parse_args()))
    except (CaseError,OSError,ValueError,KeyError) as exc: parser.exit(2,str(exc)+'\n')
    print(json.dumps({k:report[k] for k in ('complete_orbit_case','convergence_confirmed','duration_hours')}))


if __name__=='__main__': main()
