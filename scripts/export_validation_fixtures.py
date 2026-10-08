"""Export explicitly offline fixtures and planned HIL cases, never hardware data."""
import csv,hashlib,itertools,json
from pathlib import Path
import numpy as np
import plant_validation as plant
import thermal_validation as thermal
from audit_capture import TEMP_COLUMNS

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'validation'

def main():
    OUT.mkdir(exist_ok=True);av,ar=thermal.coefficients()
    for tag in ['cold_healthy','hot_aged']:
        data=np.load(ROOT/f'data/thermal_{tag}.npz')
        current=data['current'];levels=np.vstack([data['states'],data['states'][-1]])
        rows=np.column_stack([np.arange(len(current)),np.arange(len(current))*plant.CFG['Ts'],current,levels,np.full(len(current),300.),data['junction_C']])
        path=OUT/f'offline_{tag}.csv'
        with path.open('w',newline='',encoding='utf-8') as handle:
            writer=csv.writer(handle);writer.writerow(['sample_index','time_s','i_a_A','i_b_A','i_c_A','q_a','q_b','q_c','Vdc_V']+TEMP_COLUMNS);writer.writerows(rows)
        meta=dict(schema_version=1,evidence_kind='offline_simulation',source_method='generated_offline',
                  real_time_execution=False,hardware_id=None,topology='Ttype',sampling_period_s=plant.CFG['Ts'],
                  software_version='P1 offline thermal solver revision 2026-10-07',controller_provenance='Python nominal FCS current control; not external hardware',
                  calibration_provenance='Illustrative thermal coefficients in thermal_validation.py; not measured calibration',
                  machine_reference_provenance='Exact declared simulation L/E parameters',independent_machine_scale_reference=True,
                  nominal_machine={key:plant.CFG[key] for key in ['L','Rs','E','omega','td']},
                  nominal_current_noise_SD_A=plant.CFG['sigma'],actual_simulated_current_noise_SD_A=0.,reference_temperature_C=25.,
                  temperature_basis='model_junction_truth',alpha_v_V_K=av.tolist(),alpha_r_ohm_K=ar.tolist(),
                  relative_offset_drift_bound_K=0.,coefficient_relative_error_bound=0.,
                  truth_contrast=[.05,.0002] if tag=='hot_aged' else [0.,0.],
                  csv_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                  note='The last row q is ignored: each preceding endpoint row specifies its outgoing control-period level. This capture is an offline oracle fixture, not HIL or a physical temperature measurement.')
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n',encoding='utf-8')
    columns=['device_family','current_A','temperature_C','healthy_forward_drop_V','provenance']
    with (OUT/'illustrative_forward_drop_lut.csv').open('w',newline='',encoding='utf-8') as handle:
        writer=csv.writer(handle);writer.writerow(columns)
        for name,j in [('transistor',0),('diode',1)]:
            v,r=plant.nominal_devices(8)
            for temp,current in itertools.product([25,50,75,100],[0,.5,1,5,10,20,30,50]):
                drop=v[j]+av[j]*(temp-25)+current*(r[j]+ar[j]*(temp-25))
                writer.writerow([name,current,temp,drop,'illustrative affine model; not datasheet or bench calibration'])
    cases=[]
    for ambient,drift,temperature in itertools.product([25,75],['healthy','50mV_plus_0.2mohm'],['oracle_junction','case_proxy','junction_offset_drift_2K','junction_offset_drift_5K','coefficients_10percent','temperature_lag_20ms']):
        cases.append(dict(case_id=f'T{ambient}_{drift}_{temperature}',evidence_status='planned_not_executed_on_HIL',
                          ambient_C=ambient,drift=drift,temperature_intervention=temperature))
    cases.extend([dict(case_id=name,evidence_status='planned_not_executed_on_HIL',intervention=description) for name,description in [
        ('machine_mismatch','1% independent L and EMF perturbations'),('machine_reference','Independent L reference, E amplitude as nuisance'),
        ('scale_unrestricted','Both L and E unknown: test rank and reject unsupported targets'),('npc_equivalent','NPC multiple-device equivalent pair with matched temperatures'),
        ('npc_single_control','NPC distinguishable single-device pair'),('sampling_delay','Measured ADC/gate alignment and one-period delay'),
        ('solver_refinement','Halve electrical solver timestep and compare target error'),('temperature_variation','Dynamic sink variation and separate converter/machine temperature records')]])
    plan=dict(status='specification_only_no_HIL_or_bench_execution',controller_period_s=50e-6,blanking_interval_s=1e-6,
              proposed_electrical_step_s=250e-9,thermal_reporting_period_s=.001,
              note='The electrical timestep is a proposed test point, not verified device capability. Confirm platform compilation and resolve actual gate edges; refine before scoring health.',
              cases=cases,conditional_bias_budget=dict(voltage_V=.01,resistance_ohm=.00005),
              required_temperature_reference='Calibrated junction estimate with uncertainty; case temperature alone is insufficient.',
              source_urls=['https://www.typhoon-hil.com/documentation/typhoon-hil-software-manual/References/t-type_leg.html',
                           'https://www.typhoon-hil.com/documentation/typhoon-hil-software-manual/References/igbt_leg.html'])
    (OUT/'planned_hil_cases.json').write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(offline_fixtures=2,planned_HIL_cases=len(cases),hardware_runs=0),indent=2))

if __name__=='__main__':main()
