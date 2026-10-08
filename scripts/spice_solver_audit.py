"""Declared solver-tolerance intervention following full-window failures.

The physical circuit and replayed gates are unchanged. No incomplete export
is scored. The separate pair changes only voltage tolerance from the original
10 nV to 1 uV, with maximum steps of 250 and 125 ns.
"""
import json
import shutil
import numpy as np
import spice_validation as spice

def main():
    states=np.load(spice.ROOT/'data/thermal_cold_healthy.npz')['states'][:1024]
    records={};cases=[];failures=[]
    for name,step in [('aged75_tol1uV250',250e-9),('aged75_tol1uV125',125e-9)]:
        try:
            record,case=spice.run_case(states,name,75,True,step,reuse=True,vntol=1e-6)
        except RuntimeError as error:
            failed=spice.ROOT/'validation/spice_failed'
            for suffix in ['.cir','.log']:
                shutil.copy2(spice.OUT/(name+suffix),failed/(name+'_failed'+suffix))
            failures.append(dict(name=name,error=str(error),max_step_s=step,voltage_tolerance_V=1e-6))
            break
        records[name]=record;cases.append(case)
    result=dict(kind='offline_numerical_solver_intervention',hardware_run=False,
        periods=1024,voltage_tolerance_V=1e-6,cases=cases,failures=failures,
        completed_pair=len(cases)==2)
    if result['completed_pair']:
        result['max_current_difference_A']=float(abs(records[cases[1]['name']]['current']-records[cases[0]['name']]['current']).max())
        result['target_difference_V_ohm']=(np.array(cases[1]['estimate_V_ohm'])-cases[0]['estimate_V_ohm']).tolist()
        original=json.loads((spice.ROOT/'data/spice_results.json').read_text())
        baseline=next(c for c in original['cases'] if c['name']=='aged75')
        result['tolerance_only_target_difference_V_ohm']=(np.array(cases[0]['estimate_V_ohm'])-baseline['estimate_V_ohm']).tolist()
    (spice.ROOT/'data/spice_solver_audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print('SOLVER_AUDIT',json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
