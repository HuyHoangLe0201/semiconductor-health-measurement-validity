"""Generate external circuit figures and manuscript numbers from saved results."""
from pathlib import Path
import json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data';FIG=ROOT/'figures'

def main():
    result=json.loads((DATA/'spice_results.json').read_text())
    cases={x['name']:x for x in result['cases']}
    assert result['periods']==1024
    assert len(cases)==4
    refinement=result['refinement']
    assert refinement['periods']==512
    audit=json.loads((DATA/'spice_solver_audit.json').read_text())
    for case in list(cases.values())+refinement['cases']+audit['cases']:
        assert case['rank']==36 and max(case['target_defect'])<1e-8
        assert not case['log_warning_lines']
        netlist=ROOT/'validation/spice'/(case['name']+'.cir')
        assert hashlib.sha256(netlist.read_bytes()).hexdigest()==case['netlist_sha256']
        log=netlist.with_suffix('.log').read_text(errors='replace').lower()
        assert 'no. of data rows' in log and 'aborted' not in log and 'timestep too small' not in log
    for failure in audit['failures']:
        log=(ROOT/'validation/spice_failed'/(failure['name']+'_failed.log')).read_text(errors='replace').lower()
        assert 'aborted' in log and 'timestep too small' in log
    contrasts=result['contrasts']
    refined=np.array(refinement['target_difference_V_ohm'])
    result['dc_path_validation']=json.loads((ROOT/'validation/external_evidence_checks.json').read_text())['injection_path_checks']
    (DATA/'spice_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    rows=['\\begin{table}[t]\\centering\\small',
      '\\caption{Contrasts from independent nonlinear circuit replay. Both records use identical replayed gates. No measurement noise or measured device calibration is represented.}\\label{tab:spice}',
      '\\begin{tabular}{@{}lrrrr@{}}\\toprule',
      'Comparison&Injected $v$&Estimated $v$&Injected $r$&Estimated $r$\\\\',
      '&\\multicolumn{2}{c}{mV}&\\multicolumn{2}{c}{m$\\Omega$}\\\\\\midrule']
    labels=['Injected/healthy, 25 $^\\circ$C','Healthy 75/25 $^\\circ$C','Injected/healthy, 75 $^\\circ$C']
    for label,(key,values) in zip(labels,contrasts.items()):
        dv,dr=(0.,0.) if key.startswith('healthy') else (50.,.2)
        rows.append(f'{label}&{dv:.0f}&{values[0]*1000:.3f}&{dr:.1f}&{values[1]*1000:.3f}\\\\')
    rows+=['\\bottomrule\\end{tabular}','\\end{table}',
      'All four coarse records retain a rank-36 operator, with supported target-map defects below $10^{-8}$. Twelve DC operating-point probes check the six signed conduction paths in healthy and injected circuits. At a forced 10 A, only the positive upper-level path receives the expected 52 mV injected drop, agreeing within 0.001 mV. These are software circuit checks.',
      f'The maximum absolute sum of phase currents across the four full-window runs is {max(x["neutral_current_sum_max_A"] for x in cases.values())*1e6:.2f} $\\mu$A. A separate 512-period (25.6 ms) prefix of the 75 $^\\circ$C injected circuit is run at maximum internal steps of 250 and 125 ns. The maximum sampled phase-current difference is {refinement["max_current_difference_A"]*1000:.4f} mA; its fitted target changes by {abs(refined[0])*1000:.4f} mV and {abs(refined[1])*1000:.5f} m$\\Omega$. These figures concern the prefix fits, not the full-window contrasts in Table~\\ref{{tab:spice}}.',
      'Both full-window refinement attempts, at 125 and 100 ns, stopped on a transient-convergence failure at 31.4499 ms. Failed netlists and logs are retained separately; no estimates from their incomplete exports are scored. The completed prefix comparison tests local numerical sensitivity only. Full-window step convergence remains unverified, and physical model accuracy is not established by numerical agreement.']
    coarse_audit=next(x for x in audit['cases'] if x['name']=='aged75_tol1uV250')
    delta=np.array(coarse_audit['estimate_V_ohm'])-cases['aged75']['estimate_V_ohm']
    assert not audit['completed_pair'] and len(audit['failures'])==1
    rows.append(f'A separate voltage-tolerance intervention changes \\texttt{{vntol}} from $10^{{-8}}$ to $10^{{-6}}$ V while retaining the circuit, gates and other solver settings. The completed 250 ns full-window fit changes by {abs(delta[0])*1e6:.3f} $\\mu$V and {abs(delta[1])*1e6:.5f} $\\mu\\Omega$. The paired 125 ns attempt still fails at the same commutation. This tolerance check does not repair the full-window step-convergence limit.')
    (ROOT/'spice_computed_summary.tex').write_text('\n'.join(rows)+'\n',encoding='ascii')
    healthy=np.load(DATA/'spice_healthy25.npz');aged=np.load(DATA/'spice_aged25.npz')
    fig=plt.figure(figsize=(8.4,5.5),constrained_layout=True);grid=fig.add_gridspec(2,2,height_ratios=[1.3,1])
    ax=fig.add_subplot(grid[0,:]);mask=healthy['time']<=.04
    for j,phase in enumerate('abc'):ax.plot(healthy['time'][mask]*1000,healthy['current'][mask,j],lw=.9,label=f'Phase {phase}')
    ax.set(xlabel='Time (ms)',ylabel='Simulated current (A)',title=r'Independent ngspice T-type circuit: healthy 25 $^\circ$C');ax.legend(ncol=3,loc='upper right');ax.grid(alpha=.2)
    names=['Injection\n25 $^\\circ$C','Healthy\n75 vs 25 $^\\circ$C','Injection\n75 $^\\circ$C']
    for col,(unit,component,true_values) in enumerate([('Threshold contrast (mV)',0,[50,0,50]),(r'Slope contrast (m$\Omega$)',1,[.2,0,.2])]):
        ax=fig.add_subplot(grid[1,col]);values=np.array(list(contrasts.values()))[:,component]*1000
        ax.bar(np.arange(3),values,color=['#2166ac','#d95f02','#2166ac'],width=.55)
        ax.scatter(np.arange(3),true_values,marker='x',s=45,color='black',label='Injected reference',zorder=4)
        ax.set(xticks=np.arange(3),xticklabels=names,ylabel=unit);ax.axhline(0,color='.4',lw=.7);ax.grid(axis='y',alpha=.2)
    fig.savefig(FIG/'external_spice.pdf');fig.savefig(FIG/'external_spice.png',dpi=200);plt.close(fig)
    print(json.dumps({'contrasts':contrasts,'prefix_refinement':refinement},indent=2))

if __name__=='__main__':main()
