"""Offline coupled electrical/thermal study, with declared illustrative parameters.

This is not HIL, bench data, a datasheet fit, or a calibrated aging experiment.
Temperatures affect the electrical voltage drops, not just a loss postprocessor.
"""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
import plant_validation as plant

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data';FIG=ROOT/'figures'
SEED=20261010
MODEL=dict(reference_C=25.,alpha_v_transistor_V_K=-.0015,
           alpha_v_diode_V_K=-.0025,alpha_r_transistor_ohm_K=.00008,
           alpha_r_diode_ohm_K=.00005,Rj_transistor_K_W=.6,
           Rj_diode_K_W=.9,Cj_transistor_J_K=.18,Cj_diode_J_K=.12,
           Rsink_K_W=.1,Csink_J_K=20.,switching_heat_included=False)

def coefficients():
    transistor=np.arange(8)%2==0
    return (np.tile(np.where(transistor,MODEL['alpha_v_transistor_V_K'],MODEL['alpha_v_diode_V_K']),3),
            np.tile(np.where(transistor,MODEL['alpha_r_transistor_ohm_K'],MODEL['alpha_r_diode_ohm_K']),3))

def run(ambient=25.,aged=False,refine=2,N=None,zero_coefficients=False,replay=None):
    cfg=dict(plant.CFG)
    if N is not None:cfg['N']=N
    A=plant.incidence('Ttype');v0,r0=plant.nominal_devices(8)
    av,ar=coefficients()
    if zero_coefficients:av*=0;ar*=0
    av=av.reshape(3,8);ar=ar.reshape(3,8)
    dv=np.zeros((3,8));dr=dv.copy()
    if aged:dv[0,0]=.05;dr[0,0]=.0002
    tj_R=np.tile(np.where(np.arange(8)%2==0,MODEL['Rj_transistor_K_W'],MODEL['Rj_diode_K_W']),(3,1))
    tj_C=np.tile(np.where(np.arange(8)%2==0,MODEL['Cj_transistor_J_K'],MODEL['Cj_diode_J_K']),(3,1))
    dt=cfg['Ts'];count=cfg['N'];states=np.zeros((count,3),int)
    # i(3), junction T(24), common sink T, integral input heat, integral cooling.
    trajectory=np.zeros((count+1,30));trajectory[0,3:28]=ambient
    old=np.zeros(3,int)
    def derivative(time,state,levels):
        current=state[:3];tj=state[3:27].reshape(3,8);sink=state[27]
        v=v0+dv+av*(tj-25.);r=r0+dr+ar*(tj-25.)
        ids=plant.path_ids(levels,current,'Ttype');active=A[ids]
        forward=np.sum(active*(v+np.abs(current[:,None])*r),axis=1)
        pole=cfg['Vdc']/2*levels-np.sign(current)*forward
        di=(plant.P0@(pole-plant.emf(time,cfg))-cfg['Rs']*current)/cfg['L']
        power=active*(np.abs(current[:,None])*v+current[:,None]**2*r)
        flow=(tj-sink)/tj_R
        cooling=(sink-ambient)/MODEL['Rsink_K_W']
        return np.r_[di,((power-flow)/tj_C).ravel(),
                     (flow.sum()-cooling)/MODEL['Csink_J_K'],power.sum(),cooling]
    def advance(time,state,levels,h,depth=0):
        k1=derivative(time,state,levels);s2=state+h*k1/2;k2=derivative(time+h/2,s2,levels)
        s3=state+h*k2/2;k3=derivative(time+h/2,s3,levels)
        s4=state+h*k3;k4=derivative(time+h,s4,levels)
        end=state+h*(k1+2*k2+2*k3+k4)/6
        signs=np.sign(np.stack([state[:3],s2[:3],s3[:3],s4[:3],end[:3]]))
        if depth<10 and np.any(signs.min(axis=0)!=signs.max(axis=0)):
            first=advance(time,state,levels,h/2,depth+1)
            return advance(time+h/2,first,levels,h/2,depth+1)
        return end
    for k in range(count):
        time=k*dt;state=trajectory[k].copy();current=state[:3]
        if replay is None:
            candidates=plant.STATES[np.all(np.abs(plant.STATES-old)<=1,axis=1)]
            ip=np.broadcast_to(current,(len(candidates),3))
            ids=plant.path_ids(candidates,ip,'Ttype')
            drop=np.sign(ip)*(A@v0)[ids]+ip*(A@r0)[ids]
            prediction=ip+dt/cfg['L']*((cfg['Vdc']/2*candidates-drop-plant.emf(time,cfg))@plant.P0-cfg['Rs']*ip)
            cost=np.sum((prediction-plant.reference(time+dt))**2,axis=1)+.005*np.sum((candidates-old)**2,axis=1)
            levels=candidates[np.argmin(cost)]
        else:levels=replay[k]
        free=np.where(current>=0,np.minimum(old,levels),np.maximum(old,levels))
        offset=0.
        for q,duration in [(free,cfg['td']),(levels,dt-cfg['td'])]:
            steps=1 if duration==cfg['td'] else refine
            for step in range(steps):
                state=advance(time+offset,state,q,duration/steps)
                offset+=duration/steps
        state[:3]=plant.P0@state[:3];trajectory[k+1]=state
        states[k]=levels;old=levels
    energy=(trajectory[-1,3:27]-ambient)@tj_C.ravel()+MODEL['Csink_J_K']*(trajectory[-1,27]-ambient)
    energy_defect=abs(energy-(trajectory[-1,28]-trajectory[-1,29]))
    assert energy_defect<1e-8 and np.max(np.abs(trajectory[:,:3].sum(axis=1)))<1e-12
    assert np.min(trajectory[:,3:28])>=ambient-1e-8 and np.max(trajectory[:,3:28])<150
    return dict(name='Ttype',cfg=cfg,current=trajectory[:,:3],states=states,
                junction=trajectory[:,3:27],sink=trajectory[:,27],
                heat_J=trajectory[-1,28],cooling_J=trajectory[-1,29],
                energy_defect_J=float(energy_defect),ambient=ambient)

def prepare(record):
    X,y,keep=plant.observations(record)
    est,se,rank,defect,kw,details=plant.fit(X,y,keep=keep)
    assert rank==36 and np.max(defect)<1e-7
    idx=np.flatnonzero(keep);n=len(idx)
    c1=plant.CFG['L']/plant.CFG['Ts']+plant.CFG['Rs']/2
    c0=-plant.CFG['L']/plant.CFG['Ts']+plant.CFG['Rs']/2
    diagonal=(c1*c1+c0*c0)*plant.CFG['sigma']**2
    off=c1*c0*plant.CFG['sigma']**2
    diag=np.empty(n);sub=np.zeros(n);diag[0]=np.sqrt(diagonal)
    for k in range(1,n):
        sub[k]=off/diag[k-1] if idx[k]==idx[k-1]+1 else 0
        diag[k]=np.sqrt(diagonal-sub[k]**2)
    raw=kw.reshape(2,n,2).copy()
    for k in range(n-1,-1,-1):
        if k+1<n:raw[:,k]-=sub[k+1]*raw[:,k+1]
        raw[:,k]/=diag[k]
    K=np.zeros((2,len(keep),2));K[:,keep]=raw
    assert np.allclose(np.einsum('hkc,kc->h',K,y.reshape(-1,2)),est,atol=1e-9)
    shaped=X.reshape(-1,2,50)
    Gjv=np.einsum('hkc,kcj->hkj',K,shaped[:,:,:24])
    Gjr=np.einsum('hkc,kcj->hkj',K,shaped[:,:,24:48])
    return dict(X=X,y=y,keep=keep,estimate=est,se=se,rank=rank,details=details,K=K,
                Gv=Gjv,Gr=Gjr,temperature=(record['junction'][:-1]+record['junction'][1:])/2)

def corrected(design,temp=None,av=None,ar=None):
    if temp is None:return design['estimate'].copy()
    if av is None:av,ar=coefficients()
    delta=temp-25.
    return design['estimate']-np.einsum('hkj,kj,j->h',design['Gv'],delta,av)-np.einsum('hkj,kj,j->h',design['Gr'],delta,ar)

def lagged(record,tau=.02):
    true=record['junction'];measured=np.empty_like(true);measured[0]=true[0]
    pole=np.exp(-record['cfg']['Ts']/tau)
    for k in range(1,len(true)):measured[k]=pole*measured[k-1]+(1-pole)*true[k]
    # Causal 1 kHz sample-and-hold.
    sampled=measured[(np.arange(len(true))//20)*20]
    return (sampled[:-1]+sampled[1:])/2

def calibration_mc(cold,hot,offset_bound=2.,trials=1000):
    av,ar=coefficients();rng=np.random.default_rng(SEED+int(offset_bound))
    base=corrected(hot,hot['temperature'])-corrected(cold,cold['temperature'])
    Jv=np.einsum('hkj,kj->hj',hot['Gv'],hot['temperature']-25.)-np.einsum('hkj,kj->hj',cold['Gv'],cold['temperature']-25.)
    Jr=np.einsum('hkj,kj->hj',hot['Gr'],hot['temperature']-25.)-np.einsum('hkj,kj->hj',cold['Gr'],cold['temperature']-25.)
    bv=hot['Gv'].sum(axis=1);br=hot['Gr'].sum(axis=1)
    nblocks=(len(hot['temperature'])+19)//20
    def blocks(gain):
        padded=np.pad(gain,((0,0),(0,nblocks*20-gain.shape[1]),(0,0)))
        return padded.reshape(2,nblocks,20,24).sum(axis=2)
    noise_gains=[(blocks(d['Gv']),blocks(d['Gr'])) for d in [cold,hot]]
    samples=[]
    for first in range(0,trials,50):
        n=min(50,trials-first)
        cv=av*(1+rng.uniform(-.1,.1,(n,24)));cr=ar*(1+rng.uniform(-.1,.1,(n,24)))
        offset=rng.uniform(-offset_bound,offset_bound,(n,24))
        value=np.tile(base,(n,1))-(cv-av)@Jv.T-(cr-ar)@Jr.T
        value-=np.einsum('hj,nj,nj->nh',bv,cv,offset)+np.einsum('hj,nj,nj->nh',br,cr,offset)
        for sign,(gv,gr) in zip([1,-1],noise_gains):
            noise=rng.normal(0,.5,(n,nblocks,24))
            value+=sign*(np.einsum('hbj,nj,nbj->nh',gv,cv,noise)+np.einsum('hbj,nj,nbj->nh',gr,cr,noise))
        samples.append(value)
    samples=np.concatenate(samples)
    se=np.sqrt(hot['se']**2+cold['se']**2)
    return dict(trials=trials,offset_drift_bound_K=offset_bound,coefficient_relative_bound=.1,
                temperature_noise_SD_K=.5,temperature_sample_interval_s=.001,
                scope='Fixed noiseless currents and gates; exact conditional linear target map. Temperature and calibration perturbations only; coefficient errors shared between records; offsets drift only in aged record.',
                mean=samples.mean(axis=0).tolist(),SD=samples.std(axis=0,ddof=1).tolist(),
                interval95=np.quantile(samples,[.025,.975],axis=0).tolist(),
                current_noise_only_nominal_coverage=np.mean(np.abs(samples-[.05,.0002])<=1.95996398454*se,axis=0).tolist()),samples

def main():
    # Meaningful limit check: turning thermal coefficients off recovers the old plant.
    limit=run(N=128,zero_coefficients=True)
    old=plant.run(cfg=dict(plant.CFG,N=128))
    assert np.max(np.abs(limit['current']-old['current']))<1e-10
    records={};designs={}
    for tag,ambient,aged in [('cold_healthy',25.,False),('hot_healthy',75.,False),('hot_aged',75.,True)]:
        record=run(ambient,aged);records[tag]=record;designs[tag]=prepare(record)
        np.savez_compressed(DATA/f'thermal_{tag}.npz',current=record['current'],states=record['states'],junction_C=record['junction'],sink_C=record['sink'],X=designs[tag]['X'],y=designs[tag]['y'],retained=designs[tag]['keep'])
        print(tag,designs[tag]['estimate'],np.max(record['junction']),flush=True)
    cold=designs['cold_healthy'];healthy=designs['hot_healthy'];aged=designs['hot_aged']
    fine=run(75.,True,refine=4,replay=records['hot_aged']['states']);fine_design=prepare(fine)
    refinement=dict(max_current_difference_A=float(np.max(np.abs(fine['current']-records['hot_aged']['current']))),
                    max_junction_difference_K=float(np.max(np.abs(fine['junction']-records['hot_aged']['junction']))),
                    oracle_target_difference=np.abs(corrected(fine_design,fine_design['temperature'])-corrected(aged,aged['temperature'])).tolist())
    av,ar=coefficients();base_raw=corrected(cold);base_oracle=corrected(cold,cold['temperature'])
    cold_case=(records['cold_healthy']['sink'][:-1]+records['cold_healthy']['sink'][1:])/2
    hot_case=(records['hot_aged']['sink'][:-1]+records['hot_aged']['sink'][1:])/2
    drifted=aged['temperature'].copy();drifted[:,0]+=2.
    cases=[('healthy_no_compensation',healthy,None,base_raw,[0.,0.]),
           ('aged_no_compensation',aged,None,base_raw,[.05,.0002]),
           ('healthy_oracle_junction',healthy,healthy['temperature'],base_oracle,[0.,0.]),
           ('aged_oracle_junction',aged,aged['temperature'],base_oracle,[.05,.0002]),
           ('aged_sink_only',aged,np.broadcast_to(hot_case[:,None],aged['temperature'].shape),corrected(cold,np.broadcast_to(cold_case[:,None],cold['temperature'].shape)),[.05,.0002]),
           ('aged_S1_temperature_drift_2K',aged,drifted,base_oracle,[.05,.0002]),
           ('aged_coefficients_10percent',aged,aged['temperature'],corrected(cold,cold['temperature'],1.1*av,1.1*ar),[.05,.0002]),
           ('aged_temperature_lag_20ms',aged,lagged(records['hot_aged']),corrected(cold,lagged(records['cold_healthy'])),[.05,.0002])]
    rows=[]
    for name,design,temp,baseline,truth in cases:
        scaled=name=='aged_coefficients_10percent'
        estimate=corrected(design,temp,1.1*av if scaled else av,1.1*ar if scaled else ar)-baseline
        rows.append(dict(name=name,estimate=estimate.tolist(),truth=truth,bias=(estimate-truth).tolist(),nominal_SE=np.sqrt(design['se']**2+cold['se']**2).tolist()))
    calibration=[]
    for bound in [2.,5.]:
        result,samples=calibration_mc(cold,aged,bound);calibration.append(result)
        np.savez_compressed(DATA/f'thermal_calibration_{int(bound)}K.npz',estimates=samples)
    jb_v=aged['Gv'].sum(axis=1);jb_r=aged['Gr'].sum(axis=1)
    gain=jb_v*av+jb_r*ar
    l=np.zeros((2,48));l[0,0]=1;l[0,7]=-1;l[1,24]=1;l[1,30]=-1
    Gamma=np.vstack([np.diag(av),np.diag(ar)])
    assert np.allclose(gain,l@Gamma,atol=1e-10)
    budget=np.abs(gain).sum(axis=1)
    output=dict(evidence_kind='offline_simulation',seed=SEED,thermal_model=MODEL,
                thermal_solver='Fully coupled RK4 electrical and thermal states; stage-sign bisection to depth 10; heat quadrature; same FCS controller as the isothermal study.',
                scope='Illustrative conduction-heated star RC network with a common sink; no switching heat, thermal device calibration, HIL, bench, thermal interface aging, mechanical states or parasitic commutation.',
                records={tag:dict(ambient_C=rec['ambient'],junction_range_C=[float(rec['junction'].min()),float(rec['junction'].max())],sink_final_C=float(rec['sink'][-1]),input_heat_J=float(rec['heat_J']),cooling_J=float(rec['cooling_J']),energy_balance_defect_J=rec['energy_defect_J'],rank=designs[tag]['rank'],retained_periods=designs[tag]['details']['retained_periods']) for tag,rec in records.items()},
                conditions=rows,calibration_MC=calibration,integration_refinement=refinement,
                temperature_offset_gain=gain.tolist(),worst_case_offset_gain=[float(x) for x in budget],
                proposed_bias_budget=[.01,.00005],allowable_relative_offset_K=(np.array([.01,.00005])/budget).tolist(),
                zero_thermal_coefficient_current_difference_A=float(np.max(np.abs(limit['current']-old['current']))))
    print(json.dumps(output,indent=2),flush=True)
    (DATA/'thermal_results.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
    render(output,records);write_tex(output)

def render(result,records):
    t=np.arange(plant.CFG['N']+1)*plant.CFG['Ts']
    fig,axes=plt.subplots(1,2,figsize=(7.2,3),layout='constrained')
    for ax,tag in zip(axes,['cold_healthy','hot_aged']):
        rec=records[tag]
        for index,label in [(0,'S1a'),(6,'S4a'),(7,'D4a')]:ax.plot(t,rec['junction'][:,index],label=label)
        ax.plot(t,rec['sink'],'--',color='gray',label='Common sink')
        ax.set(xlabel='Time (s)',ylabel=r'Temperature ($^\circ$C)',title='Cold healthy' if tag=='cold_healthy' else 'Hot with injected drift');ax.legend(fontsize=7);ax.grid(alpha=.2)
    plant.save(fig,'thermal_trace')
    labels=['Healthy\nno T','Aged\nno T','Healthy\noracle Tj','Aged\noracle Tj','Aged\nsink only','Aged\nS1 +2 K','Aged\nslopes +10%','Aged\nT lag']
    fig,axes=plt.subplots(2,1,figsize=(7.2,4.8),layout='constrained')
    for j,ax in enumerate(axes):
        rows=result['conditions']
        ax.errorbar(np.arange(len(rows)),[1000*r['estimate'][j] for r in rows],yerr=[1.96*1000*r['nominal_SE'][j] for r in rows],fmt='o',capsize=3,label='Estimate / nominal precision')
        ax.scatter(np.arange(len(rows)),[1000*r['truth'][j] for r in rows],marker='x',color='#c05a36',label='Injected contrast')
        ax.set_xticks(range(len(labels)),labels,fontsize=7)
        ax.set_ylabel('Threshold drift (mV)' if j==0 else r'Resistance drift (m$\Omega$)');ax.legend(fontsize=7);ax.grid(alpha=.2)
    plant.save(fig,'thermal_bias')

def write_tex(result):
    labels=['Healthy, no temperature compensation','Aged, no temperature compensation','Healthy, oracle junction temperature','Aged, oracle junction temperature','Aged, common-sink proxy','Aged, $S_{1a}$ sensor drift +2 K','Aged, shared coefficient error +10\\%','Aged, 20 ms temperature lag']
    table=['\\begin{table}[t]\\centering\\small',
           '\\caption{Offline thermal interventions relative to the cold healthy record. The aged cases inject 50 mV and 0.2 m$\\Omega$; healthy cases inject zero. All coefficients and RC parameters are illustrative.}\\label{tab:thermal}',
           '\\begin{tabular}{@{}lrr@{}}\\toprule Condition&$\\widehat h_v$ (mV)&$\\widehat h_r$ (m$\\Omega$)\\\\\\midrule']
    for label,row in zip(labels,result['conditions']):table.append(label+'&'+'&'.join(f'{1000*x:.3f}' for x in row['estimate'])+'\\\\')
    table+=['\\bottomrule\\end{tabular}\\end{table}']
    rawhealthy=result['conditions'][0];oracle=result['conditions'][3];drift=result['conditions'][5]
    paragraph=f"Changing the initial sink and ambient temperature from 25 to 75 $^\\circ$C without any injected degradation gives a relative threshold estimate of {rawhealthy['estimate'][0]*1000:.3f} mV and resistance estimate of {rawhealthy['estimate'][1]*1000:.3f} m$\\Omega$ when temperature is ignored. With an oracle junction correction, the aged estimates are {oracle['estimate'][0]*1000:.3f} mV and {oracle['estimate'][1]*1000:.3f} m$\\Omega$. A 2 K drift confined to the $S_{{1a}}$ temperature channel changes these to {drift['estimate'][0]*1000:.3f} mV and {drift['estimate'][1]*1000:.3f} m$\\Omega$. The oracle is a model-truth intervention; it is not an available bench measurement. Table~\\ref{{tab:thermal}} also includes a common-sink proxy, coefficient error shared between records and a causal temperature lag."
    mc=result['calibration_MC'][1]
    refinement=result['integration_refinement']
    diagnostic=f"In {mc['trials']} conditional temperature/calibration perturbations with relative device-offset drift bounded by $\\pm5$ K, devicewise coefficient errors bounded by $\\pm10\\%$ and shared between records, and independent 0.5 K sample-and-hold noise at 1 kHz, the empirical standard deviations are {mc['SD'][0]*1000:.3f} mV and {mc['SD'][1]*1000:.3f} m$\\Omega$. Current-noise-only nominal intervals cover the injected targets in {mc['current_noise_only_nominal_coverage'][0]*100:.1f}\\% and {mc['current_noise_only_nominal_coverage'][1]*100:.1f}\\% of these perturbations. These fractions diagnose omitted calibration uncertainty; currents, gates and the electrical design are held fixed, and current noise is not regenerated. Doubling the hold subdivisions with gates replayed changes the oracle contrasts by {refinement['oracle_target_difference'][0]*1000:.4g} mV and {refinement['oracle_target_difference'][1]*1000:.4g} m$\\Omega$. The maximum thermal energy-balance defect is below $10^{{-10}}$ J. Switching heat and measured thermal parameters are absent."
    (ROOT/'thermal_computed_summary.tex').write_text('\n\n'.join(table+[paragraph,diagnostic])+'\n',encoding='utf-8')

if __name__=='__main__':main()
