"""Declared quasistatic branch-model revision, with full-window refinement.
Legacy 20-pF failures are preserved; removing Cjo changes the circuit model.
"""
import argparse,concurrent.futures,json
import numpy as np
import matplotlib.pyplot as plt
import spice_validation as s

def main(args):
    states=np.load(s.ROOT/'data/thermal_cold_healthy.npz')['states'][:1024]
    specifications=[('qs_healthy25_125',25,False,125e-9),('qs_aged25_125',25,True,125e-9),
                    ('qs_healthy75_125',75,False,125e-9),('qs_aged75',75,True,250e-9),
                    ('audit_nocap125',75,True,125e-9),('qs_aged75_62p5',75,True,62.5e-9)]
    def run(spec):return s.run_case(states,*spec,reuse=args.reuse,junction_capacitance_F=0)
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:runs=list(pool.map(run,specifications))
    records={spec[0]:pair[0] for spec,pair in zip(specifications,runs)}
    results=[pair[1] for pair in runs];lookup={r['name']:r for r in results};contrasts={}
    for key,left,right in [('age25','qs_aged25_125','qs_healthy25_125'),('age75','audit_nocap125','qs_healthy75_125'),('temperature','qs_healthy75_125','qs_healthy25_125')]:
        contrasts[key]=(np.array(lookup[left]['estimate_V_ohm'])-lookup[right]['estimate_V_ohm']).tolist()
    refinement=[]
    for coarse,fine in [('qs_aged75','audit_nocap125'),('audit_nocap125','qs_aged75_62p5')]:
        refinement.append(dict(coarse=coarse,fine=fine,periods=1024,
          max_current_difference_A=float(abs(records[fine]['current']-records[coarse]['current']).max()),
          target_difference_V_ohm=(np.array(lookup[fine]['estimate_V_ohm'])-lookup[coarse]['estimate_V_ohm']).tolist()))
    legacy=np.load(s.ROOT/'data/spice_aged75.npz')['current'];old=json.loads((s.ROOT/'data/spice_results.json').read_text())['cases']
    oldestimate=next(r['estimate_V_ohm'] for r in old if r['name']=='aged75')
    result=dict(scope='Independent offline gate replay, quasistatic Cjo=0 branch model; not HIL or switching-transient validation',
      model_revision='Only junction capacitances changed from 20 pF to zero. Gate transitions, conduction laws, temperature laws, solver tolerances and recorded gates retained.',
      periods=1024,duration_s=.0512,cases=results,contrasts=contrasts,refinement=refinement,
      legacy_comparison=dict(max_current_difference_A=float(abs(records['qs_aged75']['current']-legacy).max()),
        target_difference_V_ohm=(np.array(lookup['qs_aged75']['estimate_V_ohm'])-oldestimate).tolist()),
      full_window_completed=True,hardware_validation_performed=False)
    (s.ROOT/'data/spice_quasistatic_results.json').write_text(json.dumps(result,indent=2))
    fig,axes=plt.subplots(1,2,figsize=(7.2,3),layout='constrained')
    ax=axes[0];time=np.arange(1025)*50e-6
    for name in ['qs_aged75','audit_nocap125','qs_aged75_62p5']:
        ax.plot(time,1e6*(records[name]['current'][:,0]-records['qs_aged75_62p5']['current'][:,0]),label=f'{lookup[name]["max_step_s"]*1e9:g} ns')
    ax.set(xlabel='Time (s)',ylabel='Phase-a endpoint difference (uA)',title='Full 51.2 ms window');ax.legend(fontsize=8);ax.grid(alpha=.2)
    ax=axes[1];labels=['Age, 25 C','Age, 75 C','Healthy T change']
    for j,label in enumerate(['Threshold (mV)','Resistance (mOhm)']):ax.plot(range(3),[1000*contrasts[k][j] for k in ['age25','age75','temperature']],'o-',label=label)
    ax.set_xticks(range(3),labels,fontsize=7);ax.set(ylabel='Contrast in indicated units');ax.legend(fontsize=7);ax.grid(alpha=.2);s.plant.save(fig,'spice_quasistatic')
    tex=[r'\begin{table}[t]\centering\small',r'\caption{Quasistatic ngspice model revision ($C_{j0}=0$), with complete 1024-period records. Same-temperature injected contrasts should equal 50 mV and 0.2 m$\Omega$; the healthy temperature change has no injected degradation.}\label{tab:spiceqs}',r'\begin{tabular}{@{}lrr@{}}\toprule Comparison&$\widehat h_v$ (mV)&$\widehat h_r$ (m$\Omega$)\\\midrule']
    for label,key in zip(labels,['age25','age75','temperature']):tex.append(label+'&'+'&'.join(f'{1000*x:.3f}' for x in contrasts[key])+r'\\')
    tex+=[r'\bottomrule\end{tabular}\end{table}']
    tex.append('Full-window maximum steps of 250, 125 and 62.5 ns complete for the 75 $^\\circ$C injected circuit. The successive maximum endpoint-current differences are '+', '.join(f'{r["max_current_difference_A"]:.3g} A' for r in refinement)+'. The corresponding target changes are '+ '; '.join(f'({r["target_difference_V_ohm"][0]:.3g} V, {r["target_difference_V_ohm"][1]:.3g} $\\Omega$)' for r in refinement)+'.')
    tex.append(f'Removing the 20 pF capacitances changes the original 250 ns injected record by at most {result["legacy_comparison"]["max_current_difference_A"]:.3g} A at endpoints. This comparison measures model revision, whereas the paired step comparisons above measure numerical refinement within the revised model.')
    (s.ROOT/'spice_quasistatic_summary.tex').write_text('\n'.join(tex)+'\n')
    print('QUASISTATIC_COMPLETE',json.dumps(result),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--reuse',action='store_true');main(p.parse_args())
