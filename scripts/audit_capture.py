"""Read-only capture analysis prototype for the declared P1 T-type experiment.

Checks a neutral CSV contract; it is not a vendor driver or an authenticity
certificate. No hardware connection, energization or automatic execution.
"""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
import plant_validation as plant
import thermal_validation as thermal

DEVICES=['S1','D1','S2','D2','S3','D3','S4','D4']
TEMP_COLUMNS=[f'T_{phase}_{name}_C' for phase in 'abc' for name in DEVICES]
EVIDENCE=['offline_simulation','virtual_hil','controller_hil','bench']

def read_capture(path):
    path=Path(path);meta_path=path.with_suffix('.json')
    meta=json.loads(meta_path.read_text(encoding='utf-8'))
    if meta.get('schema_version')!=1:raise ValueError('Unsupported capture schema version.')
    if meta.get('evidence_kind') not in EVIDENCE:raise ValueError('Evidence kind must be declared explicitly.')
    if meta.get('csv_sha256')!=hashlib.sha256(path.read_bytes()).hexdigest():raise ValueError('CSV digest does not match metadata.')
    if meta.get('topology')!='Ttype':raise ValueError('This analysis prototype supports Ttype only.')
    if meta.get('evidence_kind')=='offline_simulation' and meta.get('real_time_execution'):
        raise ValueError('Offline simulation cannot be labeled real-time HIL execution.')
    if meta.get('source_method')=='generated_offline' and meta.get('evidence_kind')!='offline_simulation':
        raise ValueError('An offline-generated record cannot be relabeled as hardware evidence.')
    for field in ['software_version','controller_provenance','calibration_provenance','machine_reference_provenance']:
        if not meta.get(field):raise ValueError(f'Missing provenance: {field}.')
    if meta['evidence_kind'] in ['controller_hil','bench']:
        for field in ['hardware_id','capture_session_id','independent_reference_id']:
            if not meta.get(field):raise ValueError(f'Hardware evidence requires {field}.')
        if meta['source_method']!='hardware_capture':raise ValueError('Hardware evidence needs a hardware-capture source method.')
    with path.open(newline='',encoding='utf-8') as handle:
        reader=csv.DictReader(handle);names=reader.fieldnames
        rows=list(reader)
    required=['sample_index','time_s','i_a_A','i_b_A','i_c_A','q_a','q_b','q_c','Vdc_V']+TEMP_COLUMNS
    missing=set(required)-set(names or [])
    if missing:raise ValueError('Missing columns: '+', '.join(sorted(missing)))
    data=np.array([[float(row[name]) for name in required] for row in rows])
    if len(data)<181 or not np.isfinite(data).all():raise ValueError('Need at least 181 finite endpoint samples.')
    columns={name:data[:,k] for k,name in enumerate(required)}
    if not np.array_equal(columns['sample_index'],np.arange(len(data))):raise ValueError('Missing or duplicated sample index.')
    time=columns['time_s'];step=meta.get('sampling_period_s')
    if not np.isclose(step,plant.CFG['Ts'],rtol=0,atol=1e-12):raise ValueError('Prototype assumes the declared 50 us sampling period.')
    if abs(time[0])>1e-12 or not np.allclose(np.diff(time),step,rtol=0,atol=1e-10):raise ValueError('Timing gaps/jitter or nonzero model-phase origin; adapt the observation timing before inference.')
    q=np.column_stack([columns[f'q_{phase}'] for phase in 'abc'])
    if not np.isin(q,[-1,0,1]).all():raise ValueError('Invalid applied level.')
    if np.max(np.abs(np.diff(q[:-1],axis=0)))>1:raise ValueError('Unsupported direct +1 to -1 transition.')
    if not np.allclose(columns['Vdc_V'],plant.CFG['Vdc'],rtol=0,atol=1e-6):raise ValueError('This prototype needs adaptation for a varying/non-300 V DC link.')
    nominal=meta.get('nominal_machine',{})
    for field in ['L','Rs','E','omega','td']:
        if not np.isclose(nominal.get(field,np.nan),plant.CFG[field],rtol=0,atol=1e-12):raise ValueError('Adapt the nominal model before using a different '+field+'.')
    if not meta.get('independent_machine_scale_reference'):raise ValueError('A machine-scale reference is required for these targets.')
    if not np.isclose(meta.get('nominal_current_noise_SD_A',np.nan),plant.CFG['sigma'],rtol=0,atol=1e-12):raise ValueError('Declare the 2 mA nominal sensing reference, or adapt its covariance before inference.')
    if not np.isclose(meta.get('reference_temperature_C',np.nan),25.,rtol=0,atol=1e-12):raise ValueError('Prototype targets are defined at 25 C.')
    temps=np.column_stack([columns[name] for name in TEMP_COLUMNS])
    if not np.logical_and(temps>=0,temps<=150).all():raise ValueError('Temperature outside the declared model range.')
    if meta.get('temperature_basis') not in ['model_junction_truth','junction_estimate','measured_junction']:
        raise ValueError('Case/sink temperature cannot silently substitute for all junctions. Adapt and declare the temperature model.')
    if meta.get('temperature_basis')=='model_junction_truth' and meta['evidence_kind']=='bench':raise ValueError('Bench temperature cannot be model truth.')
    for field in ['alpha_v_V_K','alpha_r_ohm_K']:
        values=np.array(meta.get(field,[]),float)
        if values.shape!=(24,) or not np.isfinite(values).all():raise ValueError('Declare 24 finite '+field+' coefficients.')
    for field in ['relative_offset_drift_bound_K','coefficient_relative_error_bound']:
        value=meta.get(field)
        if value is None or not np.isfinite(value) or value<0:raise ValueError('Declare a nonnegative '+field+'.')
    deadline=None
    if meta['evidence_kind']=='controller_hil':
        if not meta.get('real_time_execution'):raise ValueError('Controller-HIL requires real-time execution.')
        if 'controller_latency_s' not in (names or []):raise ValueError('Controller-HIL requires measured latency per control period.')
        latency=np.array([float(row['controller_latency_s']) for row in rows[:-1]])
        if not np.isfinite(latency).all() or np.min(latency)<0:raise ValueError('Invalid measured controller latency.')
        missed=int(np.sum(latency>=step));deadline=dict(max_latency_s=float(latency.max()),missed_deadlines=missed,passed=missed==0)
        if missed:raise ValueError('Controller-HIL has missed control deadlines; do not qualify this window for health scoring.')
    record=dict(name='Ttype',cfg=dict(plant.CFG,N=len(data)-1),current=np.column_stack([columns[f'i_{p}_A'] for p in 'abc']),states=q[:-1].astype(int),junction=temps)
    common=np.max(np.abs(record['current'].sum(axis=1)))
    # Capture channels are projected to the same balanced-current basis as the paper.
    record['current']=record['current']@plant.P0
    try:design=thermal.prepare(record)
    except AssertionError as error:raise ValueError('Captured operator fails the declared rank/target-map conditions; do not attribute its unresolved targets.') from error
    av=np.array(meta['alpha_v_V_K']);ar=np.array(meta['alpha_r_ohm_K'])
    estimate=thermal.corrected(design,design['temperature'],av,ar)
    return meta,design,estimate,dict(source_csv=str(path.resolve()),evidence_kind=meta['evidence_kind'],endpoint_samples=len(data),periods=len(data)-1,rank=design['rank'],retained_periods=design['details']['retained_periods'],unprojected_current_sum_max_A=float(common),deadline=deadline)

def analyze(path,baseline=None):
    if baseline is None:raise ValueError('This prototype scores healthy-reference differences; supply --baseline.')
    meta,design,estimate,audit=read_capture(path)
    av=np.array(meta['alpha_v_V_K']);ar=np.array(meta['alpha_r_ohm_K'])
    Jv=np.einsum('hkj,kj->hj',design['Gv'],design['temperature']-25.)
    Jr=np.einsum('hkj,kj->hj',design['Gr'],design['temperature']-25.)
    se=design['se'].copy()
    if baseline:
        bmeta,bdesign,bestimate,baudit=read_capture(baseline)
        if bmeta['evidence_kind']!=meta['evidence_kind']:raise ValueError('Do not mix simulated and hardware baseline evidence.')
        if bmeta['calibration_provenance']!=meta['calibration_provenance']:raise ValueError('Shared-coefficient bound requires the same calibration provenance.')
        if not np.array_equal(bmeta['alpha_v_V_K'],meta['alpha_v_V_K']) or not np.array_equal(bmeta['alpha_r_ohm_K'],meta['alpha_r_ohm_K']):raise ValueError('Baseline and target need matching calibration coefficients.')
        estimate-=bestimate;se=np.sqrt(se**2+bdesign['se']**2)
        Jv-=np.einsum('hkj,kj->hj',bdesign['Gv'],bdesign['temperature']-25.)
        Jr-=np.einsum('hkj,kj->hj',bdesign['Gr'],bdesign['temperature']-25.)
        audit['baseline']=baudit
    gain=design['Gv'].sum(axis=1)*av+design['Gr'].sum(axis=1)*ar
    fraction=meta['coefficient_relative_error_bound']
    offset_bound=np.abs(gain).sum(axis=1)*meta['relative_offset_drift_bound_K']
    # Bound the additional coefficient/offset product instead of dropping it.
    cross_bound=fraction*(np.abs(design['Gv'].sum(axis=1))*np.abs(av)+np.abs(design['Gr'].sum(axis=1))*np.abs(ar)).sum(axis=1)*meta['relative_offset_drift_bound_K']
    slope_bound=fraction*(np.abs(Jv)*np.abs(av)+np.abs(Jr)*np.abs(ar)).sum(axis=1)
    component_bound=offset_bound+cross_bound+slope_bound
    budget=np.array([.01,.00005])
    audit.update(dict(status='offline_capture_analysis' if meta['evidence_kind']=='offline_simulation' else 'capture_analysis_requires_physical_reference_review',
                      hardware_validation_performed=False,
                      target_definition=['v_S1a_minus_v_D4a_at_25C','r_S1a_minus_r_S4a_at_25C'],
                      healthy_reference_subtracted=bool(baseline),estimate=estimate.tolist(),
                      current_noise_only_nominal_SE=se.tolist(),
                      calibration_component_bound=component_bound.tolist(),
                      calibration_component_bias_budget=budget.tolist(),
                      calibration_component_budget_met=(component_bound<=budget).tolist(),
                      limitations='Conditional prototype with declared affine model. Bounds cover constant offset drift and shared coefficient uncertainty only; dynamic temperature error, noise, commutation/model discrepancy and hardware-reference uncertainty are omitted. Provenance fields are consistency checks, not authentication. No HIL/bench execution is performed by this script.'))
    return audit

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('capture',type=Path);parser.add_argument('--baseline',type=Path,required=True);parser.add_argument('--output',type=Path)
    args=parser.parse_args()
    try:result=analyze(args.capture,args.baseline)
    except (ValueError,KeyError,FileNotFoundError) as error:
        parser.exit(2,'Capture rejected: '+str(error)+'\n')
    encoded=json.dumps(result,indent=2)
    if args.output:args.output.write_text(encoded+'\n',encoding='utf-8')
    print(encoded)

if __name__=='__main__':main()
