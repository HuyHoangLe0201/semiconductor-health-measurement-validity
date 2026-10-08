"""Independent DC conduction-path probes and external-data consistency checks."""
from pathlib import Path
import json,subprocess,re,hashlib
import numpy as np
from scipy.io import loadmat
import spice_validation as spice

ROOT=spice.ROOT;OUT=spice.OUT

def probe(q,current,aged):
    name=f'path_q{q}_i{int(current)}_'+('aged' if aged else 'healthy')
    lines=[name,'Vpos pos 0 150','Vneg neg 0 -150',
       '.model td D(Is=1e-15 N=1 Rs=0.015 Cjo=20p)',
       '.model ad D(Is=1e-13 N=1 Rs=0.015 Cjo=20p)',
       '.temp 25','.options rshunt=1e8 reltol=1e-8 abstol=1e-10',f'Iload p 0 {current}']
    for j,(collector,emitter) in enumerate([('pos','p'),('0','mid'),('p','mid'),('p','neg')],1):
        G=100 if spice.gates(q)[j-1] else 1e-6
        lines.append(f'Rsw{j} {collector} y{j} {1/G:.12g}')
        if aged and j==1:lines += ['Vshift y1 xshift 0.05','Rwire1 xshift x1 0.0102']
        else:lines += [f'Rwire{j} y{j} x{j} 0.01']
        lines += [f'Ds{j} x{j} {emitter} td',f'Da{j} {emitter} {collector} ad']
    lines += ['.control','set numdgt=15','op','print v(p) '+ ' '.join(f'@ds{j}[id] @da{j}[id]' for j in range(1,5)),'quit','.endc','.end']
    p=OUT/(name+'.cir');p.write_text('\n'.join(lines)+'\n',encoding='ascii')
    result=subprocess.run([str(spice.EXE),'-b','-o',name+'.log',p.name],cwd=OUT,capture_output=True,text=True,timeout=60)
    log=(OUT/(name+'.log')).read_text(errors='replace')
    if result.returncode:raise RuntimeError(log)
    values={m[0]:float(m[1]) for m in re.findall(r'([\w@()\[\]]+)\s*=\s*([+-]?[0-9.]+e[+-][0-9]+)',log,re.I)}
    assert 'v(p)' in values,log
    currents={f'S{j}':values[f'@ds{j}[id]'] for j in range(1,5)}
    currents.update({f'D{j}':values[f'@da{j}[id]'] for j in range(1,5)})
    expected={(1,1):['S1'],(0,1):['S2','D3'],(-1,1):['D4'],(1,-1):['D1'],(0,-1):['S3','D2'],(-1,-1):['S4']}[(q,int(np.sign(current)))]
    active=sorted(k for k,v in currents.items() if v>1.)
    assert active==sorted(expected),(name,active,expected)
    return dict(name=name,q=q,forced_current_A=current,aged=aged,pole_voltage_V=values['v(p)'],active=active,expected=sorted(expected),currents_A=currents)

def main():
    probes=[];shift_checks=[]
    for q in [1,0,-1]:
        for current in [10.,-10.]:
            h=probe(q,current,False);a=probe(q,current,True);probes.extend([h,a])
            difference=a['pole_voltage_V']-h['pole_voltage_V']
            expected=-.052 if q==1 and current>0 else 0.
            assert abs(difference-expected)<1e-4,(q,current,difference,expected)
            shift_checks.append(dict(q=q,current_A=current,delta_pole_V=difference,expected_V=expected))
    manifest=json.loads((ROOT/'external/nasa/provenance.json').read_text())
    z=np.load(ROOT/'data/nasa_onstate_blocks.npz');summary=z['summary'];blocks=z['blocks']
    assert len(summary)==310 and len(blocks)==107681
    assert 423-len(summary)==len(manifest['excluded'])
    sampling=[]
    for file in manifest['files']:
        p=ROOT/'external/nasa'/file['file']
        h=hashlib.sha256(p.read_bytes()).hexdigest();assert h==file['sha256']
        a=summary[summary[:,0]==file['device']]
        assert np.all(np.diff(a[:,2])>=0),'source captures are not chronological'
        assert np.all(a[:,12]<=2.)
        raw=loadmat(p,simplify_cells=True)['measurement']
        sampling.append(dict(device=file['device'],source_captures=len(raw['transient']),
            sample_counts=sorted({len(x['timeDomain']['collectorEmitterVoltage']) for x in raw['transient']}),
            sample_intervals_s=sorted({float(x['timeDomain']['dt']) for x in raw['transient']})))
        del raw
    # Reconstruct one source waveform through a vectorized independent blocking
    # route and compare saved medians, rather than merely re-reading CSV numbers.
    m=loadmat(ROOT/'external/nasa/Device2__1.mat',simplify_cells=True)['measurement']
    td=m['transient'][0]['timeDomain'];dt=td['dt'];width=round(2e-6/dt);guard=round(10e-6/dt)
    gate=td['gateEmitterVoltage'];i=td['collectorEmitterCurrentSignal'];v=td['collectorEmitterVoltage']
    good=(gate>10)&(i>1)&(v>0)&np.isfinite(i+v);good[:guard]=False;good[-guard:]=False
    for e in np.flatnonzero(np.diff((gate>10).astype(int)))+1:good[max(0,e-guard):min(len(i),e+guard)]=False
    mask=good[:len(good)//width*width].reshape(-1,width).all(axis=1)
    im=np.median(i[:len(mask)*width].reshape(-1,width),axis=1)[mask]
    vm=np.median(v[:len(mask)*width].reshape(-1,width),axis=1)[mask]
    saved=blocks[(blocks[:,0]==2)&(blocks[:,1]==1)]
    np.testing.assert_array_equal(im,saved[:,3]);np.testing.assert_array_equal(vm,saved[:,4])
    result=dict(evidence_kind='software_and_third_party_data_checks',hardware_run=False,
                source_hashes_match=True,counts_reconcile=True,temperature_time_match_within_2s=True,
                capture_order_chronological=True,independent_source_median_check=True,
                signed_path_checks=probes,injection_path_checks=shift_checks,raw_sampling=sampling)
    (ROOT/'validation/external_evidence_checks.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({'checks':'passed','signed_path_cases':len(probes),'injection_path_cases':len(shift_checks),'NASA_captures':len(summary),'NASA_blocks':len(blocks),'raw_sampling':sampling},indent=2),flush=True)

if __name__=='__main__':main()
