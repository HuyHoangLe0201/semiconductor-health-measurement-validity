"""Revision figures and text generated from stored, computed evidence."""
from pathlib import Path
import json
import numpy as np
import sympy as sp
from scipy.stats import beta
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import run_experiments as base

ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'data';FIG=ROOT/'figures'

def save(fig,name):
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)

def exact_heterogeneous():
    records=[(0,u,n) for u in [1,2] for n in [0,1]]+[(1,2,0),(1,2,1),(1,4,1),(1,4,2)]
    rows=[]
    for p,u,n in records:
        e=np.eye(2,dtype=int)[p];rows.append(np.r_[e,u*e,n])
    X=sp.Matrix(np.array(rows).tolist());Z=X[:,:2];H=X[:,2:]
    F=Z.T*(sp.eye(8)-H*(H.T*H).inv()*H.T)*Z
    D=sp.diag(sp.Rational(1,20),sp.Rational(1,20));b=sp.Matrix([sp.Rational(1,40),-sp.Rational(1,40)])
    k=sp.Rational(1,2)*(sp.Rational(1,2)-sp.Rational(3,4)**2/sp.Rational(5,2))+sp.Rational(1,2)*(sp.Rational(3,2)-sp.Rational(7,2)**2/10)
    derived=8*(D-b*b.T/k)
    assert F==derived
    contrast=sp.Matrix([1,-1]);level=sp.Matrix([1,1])
    values={'F':[[str(x) for x in row] for row in F.tolist()], 'penalty_direction':[str(x) for x in D.inv()*b],'contrast_bound':str((contrast.T*F.inv()*contrast)[0]),'sum_bound':str((level.T*F.inv()*level)[0]),'known_dt_contrast_bound':str((contrast.T*(8*D).inv()*contrast)[0])}
    return values

def scaling_fixture():
    data=np.load(DATA/'Ttype_design.npz');X=data['X'];A=data['incidence'];ids=data['path_ids'];i=data['currents_A']
    levels=np.array([1,0,-1,1,0,-1]);signs=np.array([1,1,1,-1,-1,-1])
    kv=np.linalg.lstsq(A,signs*levels*150,rcond=None)[0]
    v0,r0=np.where(np.arange(8)%2==0,1.1,.9),np.where(np.arange(8)%2==0,.02,.015)
    v=np.tile(v0,3);v[0]+=.05;r=np.tile(r0,3);r[0]+=.0002
    complete=np.r_[v,r,.4,.25]
    C=base.np.array([[1,-1,0],[1,1,-2]],float);C/=np.linalg.norm(C,axis=1)[:,None]
    f=np.array([np.sin(.03*k+np.array([0,-2*np.pi/3,2*np.pi/3])) for k in range(len(i))])
    cf=(f@C.T).ravel();cu=(150*levels[ids]@C.T).ravel()
    derivative=(cu+X@complete-20*cf)/.002
    extended=np.column_stack([X,-derivative,-cf])
    tangent=np.r_[v-np.tile(kv,3),r,.4,.25,.002,20]
    error=float(np.max(np.abs(extended@tangent)))
    assert error<1e-10 and base.rank(extended)==37
    np.savez_compressed(DATA/'machine_scale_fixture.npz',X=extended,tangent=tangent)
    return {'rank':base.rank(extended),'parameters':52,'tangent_residual_max_V':error,'threshold_contrast_tangent_V':float(tangent[0]-tangent[7]),'resistance_contrast_tangent_ohm':float(tangent[24]-tangent[30]),'scope':'Equivalent affine observation model with ideal current derivatives, not switching-event plant.'}

def main():
    # Recreate the three algebraic figures, then improve ambiguity visualization.
    base.ambiguity_figures()
    Xn=np.load(DATA/'NPC_design.npz')['X'];Xt=np.load(DATA/'Ttype_design.npz')['X']
    a=np.zeros(62);a[0]=a[8]=.05;b=np.zeros(62);b[2]=.05;s=a.copy();s[8]=0
    t1=np.zeros(50);t1[0]=.05;t2=np.zeros(50);t2[7]=.05
    fig,axes=plt.subplots(2,3,figsize=(7.8,3.7),layout='constrained')
    for j,(X,u,v,title) in enumerate([(Xn,a,b,'NPC: multiple devices'),(Xn,s,b,'NPC: single device'),(Xt,t1,t2,'T-type: outer devices')]):
        y1=1000*(X@u)[:100];y2=1000*(X@v)[:100]
        axes[0,j].plot(y1,label='A');axes[0,j].plot(y2,'--',label='B');axes[0,j].set_title(title);axes[0,j].legend(fontsize=7)
        axes[1,j].plot(y1-y2,color='#438565');axes[1,j].axhline(0,color='gray',lw=.6);axes[1,j].set_xlabel('Residual coordinate')
    axes[0,0].set_ylabel('Residual (mV)');axes[1,0].set_ylabel('A minus B (mV)');save(fig,'ambiguity')
    rows=json.loads((DATA/'manuscript_results.json').read_text())['physical_MC']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
    for j,ax in enumerate(axes):
        for tag,label,offset,color in [('coverage','GLS / MA(1)',-.025,'#22577a'),('naive_white_coverage','OLS / white',.025,'#c05a36')]:
            cov=np.array([row[tag][j] for row in rows]);trials=np.array([row['trials'] for row in rows]);count=np.rint(cov*trials).astype(int)
            lo=beta.ppf(.025,count,trials-count+1);hi=beta.ppf(.975,count+1,trials-count)
            ax.errorbar(np.array([row['periods'] for row in rows])*(1+offset),cov*100,yerr=np.vstack([cov-lo,hi-cov])*100,fmt='o',capsize=3,color=color,label=label)
        ax.axhline(95,color='gray',lw=.8);ax.set(xlabel='Observation periods',ylabel='Coverage (%)',title='Threshold contrast' if j==0 else 'Resistance contrast');ax.legend(fontsize=7);ax.grid(alpha=.2)
    save(fig,'coverage')
    out={'heterogeneous_fixture':exact_heterogeneous(),'machine_scale_fixture':scaling_fixture()}
    if (DATA/'plant_results.json').exists():
        plant=json.loads((DATA/'plant_results.json').read_text())
        macros=[]
        for j,tag in [(0,'Voltage'),(1,'Resistance')]:
            vals=plant['detection'][j]
            for suffix,key,scale,precision in [('Power','power',100,2),('Half','halfwidth95',1000,3),('MDE','approx_MDE80',1000,3)]:
                macros.append('\\newcommand{\\Detection'+tag+suffix+'}{'+f'{scale*vals[key]:.{precision}f}'+'}')
        macros.append('\\newcommand{\\NPCAmbiguousDifference}{'+f"{plant['NPC_gate_replay']['ambiguous_max_current_difference_A']:.2g}"+'}')
        macros.append('\\newcommand{\\NPCSingleDifference}{'+f"{plant['NPC_gate_replay']['single_device_max_current_difference_A']:.5f}"+'}')
        (ROOT/'revision_values.tex').write_text('\n'.join(macros)+'\n',encoding='utf-8')
        summary=['\\begin{table}[t]','\\caption{Switching-event plant diagnostics for the T-type health contrasts. Standard errors use the nominal measured-design Gaussian model; they exclude model discrepancy.}\\label{tab:plant}','\\centering\\small','\\begin{tabular}{@{}lrrrr@{}}\\toprule','Condition&$\\widehat h_v$ (mV)&$s_v$ (mV)&$\\widehat h_r$ (m$\\Omega$)&$s_r$ (m$\\Omega$)\\\\\\midrule']
        labels={'matched':'Matched, noiseless','noisy_feedback':'Noisy feedback and sensing','L_and_EMF_mismatch':'1\\% $L$ and $E$ mismatch','extended_nuisance':'Both $L$ and $E$ unknown','independent_L_reference':'Known $L$, unknown $E$'}
        for rec in plant['conditions']:
            if rec['name']=='extended_nuisance':
                summary.append(labels[rec['name']]+'&Rejected&--&Rejected&--\\\\')
                continue
            values=[rec['estimate'][0],rec['nominal_SE'][0],rec['estimate'][1],rec['nominal_SE'][1]]
            summary.append(labels[rec['name']]+'&'+'&'.join(f'{1000*x:.3f}' for x in values)+'\\\\')
        stationary,changed=plant['healthy_reference']['records']
        values=[stationary['relative_estimate'][0],stationary['nominal_relative_SE'][0],stationary['relative_estimate'][1],stationary['nominal_relative_SE'][1]]
        summary.append('Stationary healthy reference&'+'&'.join(f'{1000*x:.3f}' for x in values)+'\\\\')
        summary+=['\\bottomrule\\end{tabular}','\\end{table}']
        mc=plant['measurement_MC']
        summary.append(f"Table~\\ref{{tab:plant}} distinguishes finite diagnostic fits from the rejected unrestricted fit. Across {mc['trials']} independently perturbed sensor records with gates held fixed, the empirical standard deviations are {mc['SD'][0]*1000:.3f} mV and {mc['SD'][1]*1000:.3f} m$\\Omega$; the root-mean-square errors are {mc['RMSE'][0]*1000:.3f} mV and {mc['RMSE'][1]*1000:.3f} m$\\Omega$. Each nominal interval uses the standard error from its refitted measured design; the mean standard errors are {mc['mean_refitted_SE'][0]*1000:.3f} mV and {mc['mean_refitted_SE'][1]*1000:.3f} m$\\Omega$. Coverage is {mc['nominal_coverage'][0]*100:.1f}\\% and {mc['nominal_coverage'][1]*100:.1f}\\%. In particular, the threshold interval undercovers despite a small standard error. With only {mc['trials']} trials, these coverage estimates have substantial binomial uncertainty and are diagnostic rather than a calibration claim.")
        bad=next(rec for rec in plant['conditions'] if rec['name']=='extended_nuisance')
        summary.append(f"The unrestricted two-machine-nuisance fit is rejected: it drives the inferred total inductance and back-EMF amplitude essentially to zero. Its numerical rank is {bad['rank']}, compared with 37 in the exact equivalent-model fixture, and its smallest scaled retained singular value is {bad['minimum_scaled_singular_value']:.3g}. The nominal threshold standard error would be {bad['nominal_SE'][0]:.1f} V. This is numerical lifting of a scale ambiguity into a nearly singular direction, not restored health attribution. Physical constraints can exclude the zero-scale boundary but require a separate constrained inference analysis.")
        assert plant['integration_target_difference'][0]<1e-9 and plant['integration_target_difference'][1]<1e-12
        summary.append(f"The recorded base operator has rank {plant['conditions'][0]['rank']} and retains {plant['conditions'][0]['retained_periods']} periods. Replaying the base gates with doubled hold-interval integration subdivisions changes current by at most {plant['integration_max_difference_A']*1e6:.2f} $\\mu$A; the fitted contrasts change by less than $10^{{-6}}$ mV and $10^{{-9}}$ m$\\Omega$. This checks the numerical integrator; it does not bound the equivalent observation model's discrepancy. The noisy-feedback condition reruns the controller for one sensor realization. The {mc['trials']}-trial study instead perturbs measurements on a recorded trajectory and reconstructs and refits its noisy regressors each time; it is not a closed-loop Monte Carlo experiment.")
        summary.append('\\subsection{Healthy-reference intervention}\\label{sec:healthyref}')
        summary.append(f"A healthy-reference difference is a distinct observation task. Two independent healthy trajectories are generated, one with the same 1\\% machine mismatch as the aged record and one with the nominal machine before its parameters change. Subtracting the fitted healthy contrasts from the aged contrasts gives {stationary['relative_estimate'][0]*1000:.3f} mV and {stationary['relative_estimate'][1]*1000:.3f} m$\\Omega$ when the machine mismatch remains stationary. The nominal independent-record standard errors are {stationary['nominal_relative_SE'][0]*1000:.3f} mV and {stationary['nominal_relative_SE'][1]*1000:.3f} m$\\Omega$. When the machine parameters change between the healthy and aged records, the threshold estimate instead becomes {changed['relative_estimate'][0]:.3f} V, while the resistance estimate is {changed['relative_estimate'][1]*1000:.3f} m$\\Omega$. Thus a stationary baseline can cancel a large common bias, but subtracting a baseline does not certify that the remaining change is semiconductor degradation. This is an intervention within the present estimator, not a reproduced benchmark of Ou et al.'s grid-connected method \\cite{{ou2025}}.")
        (ROOT/'plant_computed_summary.tex').write_text('\n\n'.join(summary)+'\n',encoding='utf-8')
    (DATA/'revision_checks.json').write_text(json.dumps(out,indent=2),encoding='utf-8');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
