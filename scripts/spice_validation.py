"""Independent ngspice physical T-type circuit with generic conduction laws.

Offline gate replay, not controller HIL or a manufacturer-calibrated IGBT.
The circuit solver determines conduction through switches/diodes and an
isolated star neutral; it does not use the estimator's path incidence matrix.
"""
from pathlib import Path
import argparse,json,subprocess,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import plant_validation as plant

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'validation'/'spice';OUT.mkdir(parents=True,exist_ok=True)
EXE=ROOT/'external/ngspice47/Spice64/bin/ngspice_con.exe'

def gates(q):
    # P=(S1,S2), O=(S2,S3), N=(S3,S4). Retain unchanged gates
    # during 1-us blanking; defer only newly activated devices.
    return np.array([q==1,q>=0,q<=0,q==-1],int)

def pwl_gates(states,phase,Ts,td):
    old=gates(0);points=[[(0.,int(x))] for x in old]
    edge=50e-9
    for k,q in enumerate(states[:,phase]):
        new=gates(q);t=k*Ts
        for j in range(4):
            if new[j]!=old[j]:
                switch=t+(td if new[j] else 0.)
                if switch>0:points[j].append((max(0,switch-edge),int(old[j])))
                points[j].append((switch+edge,int(new[j])))
        old=new
    for j in range(4):points[j].append((len(states)*Ts,int(old[j])))
    return points

def write_netlist(states,name,temperature=25.,aged=False,maxstep=250e-9,vntol=1e-8,gate_smoothing=False,junction_capacitance_F=20e-12,branch_model='series_diode',gate_conductance='logarithmic'):
    cfg=plant.CFG;N=len(states);Ts=cfg['Ts'];td=cfg['td'];path=OUT/(name+'.cir')
    lines=[f'P1 generic nonlinear T-type circuit: {name}',
           '* Generic diode-plus-gated-switch surrogate; no device calibration.',
           'Vpos pos 0 150','Vneg neg 0 -150',
           '.model gateSW SW(Ron=0.001 Roff=1e9 Vt=0.5 Vh=0)',
           '.model transistorD D(Is=1e-15 N=1 Rs=0.015 Cjo=20p Tt=0)',
           '.model antiparallelD D(Is=1e-13 N=1 Rs=0.015 Cjo=20p Tt=0)',
           f'.temp {temperature}',
           f'.options reltol=1e-6 abstol=1e-9 vntol={vntol:.12g} method=gear maxord=2 rshunt=1e8 itl4=200',
           'Rneutral neutral 0 1e9','.save i(La) i(Lb) i(Lc)']
    for p,phase in enumerate('abc'):
        wave=pwl_gates(states,p,Ts,td)
        for j,points in enumerate(wave,1):
            lines.append(f'Vg{phase}{j} g{phase}{j} 0 PWL(')
            for start in range(0,len(points),6):
                lines.append('+ '+' '.join(f'{t:.12g} {v}' for t,v in points[start:start+6]))
            lines.append('+ )')
        # Antiparallel packages: high bus -> pole; middle -> common emitter;
        # pole -> common emitter; pole -> low bus.
        for j,(collector,emitter) in enumerate([('pos',f'p{phase}'),('0',f'm{phase}'),(f'p{phase}',f'm{phase}'),(f'p{phase}','neg')],1):
            anode=f'x{phase}{j}';switch_end=f'y{phase}{j}'
            # Smooth logarithmic conductance through a 100-ns gate transition.
            g=f'V(g{phase}{j})'
            if gate_smoothing:g=f'(3*{g}^2-2*{g}^3)'
            if branch_model=='direct_rectifier':
                vth=.9-.0015*(temperature-25)+(0.05 if aged and phase=='a' and j==1 else 0)
                resistance=.025+.00008*(temperature-25)+(0.0002 if aged and phase=='a' and j==1 else 0)
                x=f'(V({collector},{emitter})-({vth:.12g}))'
                lines += [f'Bt{phase}{j} {collector} {emitter} I=0.5*({x}+sqrt({x}^2+0.01^2))/({resistance:.12g}+1e6*exp(-18.42068074395*V(g{phase}{j})))',
                    f'Cce{phase}{j} {collector} {emitter} 20p',f'Da{phase}{j} {emitter} {collector} antiparallelD']
                continue
            law=f'1e-6*exp(18.42068074395*{g})'
            if gate_conductance=='smooth_linear':
                gate=f'V(g{phase}{j})';law=f'(1e-6+99.999999*(3*{gate}^2-2*{gate}^3))'
            lines.append(f'Bsw{phase}{j} {collector} {switch_end} I=V({collector},{switch_end})*{law}')
            if aged and phase=='a' and j==1:
                lines += [f'Bdrift {switch_end} zshift V=0.05',f'Rdrift zshift {anode} 0.0102']
            else:lines.append(f'Rwire{phase}{j} {switch_end} {anode} 0.01')
            lines += [f'Ds{phase}{j} {anode} {emitter} transistorD',f'Da{phase}{j} {emitter} {collector} antiparallelD']
        phase_deg=-120*p
        lines += [f'L{phase} p{phase} r{phase} 0.002 IC=0',f'R{phase} r{phase} e{phase} 0.25',f'Bemf{phase} e{phase} neutral V=20*sin(2*pi*25*time+({phase_deg})*pi/180)*min(time/1e-6,1)']
    # Uniform output at controller endpoints; internal step is limited separately.
    lines += ['.control','set klu','set noaskquit','set numdgt=15','set wr_singlescale','set wr_vecnames',
              f'tran {Ts:.12g} {N*Ts:.12g} 0 {maxstep:.12g}',
              'linearize i(La) i(Lb) i(Lc)',
              f'wrdata {name}_current.txt i(La) i(Lb) i(Lc)',
              'quit','.endc','.end']
    if junction_capacitance_F!=20e-12:lines=[line.replace('Cjo=20p',f'Cjo={junction_capacitance_F:.12g}') for line in lines]
    path.write_text('\n'.join(lines)+'\n',encoding='ascii')
    return path

def run_case(states,name,temperature,aged,maxstep,reuse=False,vntol=1e-8,gate_smoothing=False,junction_capacitance_F=20e-12,branch_model='series_diode',gate_conductance='logarithmic'):
    previous_path=OUT/(name+'.cir');previous=previous_path.read_bytes() if previous_path.exists() else None
    path=write_netlist(states,name,temperature,aged,maxstep,vntol=vntol,gate_smoothing=gate_smoothing,junction_capacitance_F=junction_capacitance_F,branch_model=branch_model,gate_conductance=gate_conductance)
    output=OUT/(name+'_current.txt');logpath=OUT/(name+'.log')
    log=logpath.read_text(errors='replace') if logpath.exists() else ''
    cached=reuse and previous==path.read_bytes() and output.exists() and log.rstrip().endswith('ngspice-47 done') and 'aborted' not in log.lower() and 'timestep too small' not in log.lower()
    if not cached:
        completed=subprocess.run([str(EXE),'-b','-o',name+'.log',path.name],cwd=OUT,capture_output=True,text=True,timeout=1200)
        log=logpath.read_text(errors='replace')
        if completed.returncode or not output.exists() or not log.rstrip().endswith('ngspice-47 done') or 'aborted' in log.lower() or 'timestep too small' in log.lower():raise RuntimeError(log[-6000:]+completed.stderr)
    a=np.loadtxt(output,skiprows=1)
    assert a.shape==(len(states)+1,4),a.shape
    np.testing.assert_allclose(a[:,0],np.arange(len(a))*plant.CFG['Ts'],atol=1e-11,rtol=0)
    record=dict(current=a[:,1:],states=states,name='Ttype',cfg={**plant.CFG,'N':len(states)})
    X,y,keep=plant.observations(record);estimate,se,rank,defect,_,details=plant.fit(X,y,keep=keep)
    result=dict(name=name,temperature_C=temperature,aged=aged,max_step_s=maxstep,voltage_tolerance_V=vntol,gate_smoothing=gate_smoothing,junction_capacitance_F=junction_capacitance_F,branch_model=branch_model,gate_conductance=gate_conductance,estimate_V_ohm=estimate.tolist(),
                rank=rank,target_defect=defect.tolist(),retained=details['retained_periods'],
                neutral_current_sum_max_A=float(abs(record['current'].sum(axis=1)).max()),
                netlist_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),log=path.with_suffix('.log').name,
                log_warning_lines=[l for l in log.splitlines() if 'warning' in l.lower() or 'failed' in l.lower()])
    np.savez_compressed(ROOT/'data'/f'spice_{name}.npz',time=a[:,0],current=a[:,1:],states=states)
    print('SPICE_CASE',json.dumps(result),flush=True)
    return record,result

def main(args):
    N=args.periods
    # Gates recorded by the independent existing controller, replayed unchanged.
    z=np.load(ROOT/'data/thermal_cold_healthy.npz')
    print('TRACE_KEYS',z.files,flush=True)
    states=z['states'][:N]
    version=subprocess.run([str(EXE),'--version'],capture_output=True,text=True).stdout
    results=[];records={}
    cases=[('healthy25',25,False),('aged25',25,True),('healthy75',75,False),('aged75',75,True)]
    if args.smoke:cases=cases[:1]
    for name,temp,aged in cases:
        rec,res=run_case(states,name,temp,aged,250e-9,reuse=args.reuse);records[name]=rec;results.append(res)
    if args.smoke:return
    # Full-window 125/100-ns attempts failed at 31.4499 ms. A declared prefix
    # comparison is scored separately; incomplete exports are never analysed.
    prefix=states[:512]
    refinement=[];prefix_records={}
    for name,step in [('aged75_prefix250',250e-9),('aged75_prefix125',125e-9)]:
        rec,res=run_case(prefix,name,75,True,step,reuse=args.reuse)
        prefix_records[name]=rec;refinement.append(res)
    contrasts={}
    for key,left,right in [('aged25_vs_healthy25','aged25','healthy25'),('healthy75_vs_healthy25','healthy75','healthy25'),('aged75_vs_healthy75','aged75','healthy75')]:
        l=next(r for r in results if r['name']==left);r=next(r for r in results if r['name']==right)
        contrasts[key]=(np.array(l['estimate_V_ohm'])-r['estimate_V_ohm']).tolist()
    diff=prefix_records['aged75_prefix125']['current']-prefix_records['aged75_prefix250']['current']
    target_diff=np.array(refinement[1]['estimate_V_ohm'])-refinement[0]['estimate_V_ohm']
    final=dict(evidence_kind='offline_external_circuit_simulation',hardware_validation_performed=False,
               simulator_version=version,periods=N,cases=results,contrasts=contrasts,
               refinement=dict(periods=len(prefix),duration_s=len(prefix)*plant.CFG['Ts'],cases=refinement,
                   max_current_difference_A=float(abs(diff).max()),target_difference_V_ohm=target_diff.tolist(),
                   scope='Prefix only; full-window refinement failed at 31.4499 ms for both 125 ns and 100 ns.',
                   full_window_failed_max_steps_ns=[125,100]),
               limitations='Generic nonlinear conduction surrogate; gate replay, no closed-loop external controller; no calibrated IGBT, switching parasitics or dynamic thermal network.')
    (ROOT/'data/spice_results.json').write_text(json.dumps(final,indent=2),encoding='utf-8')
    print('SPICE_RESULT',json.dumps(final,indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--periods',type=int,default=1024);p.add_argument('--smoke',action='store_true');p.add_argument('--reuse',action='store_true');main(p.parse_args())
