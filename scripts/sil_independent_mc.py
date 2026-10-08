"""Independent-record closed-loop Monte Carlo, with held-out null calibration.

All execution is offline. The original circuit, controller and estimator are
retained. The additional native harness only batches acquisition and exposes
the nominal QR target covariance. A calibrated null threshold is not a joint
likelihood or a calibrated physical uncertainty statement.
"""
from pathlib import Path
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse, ctypes as ct, hashlib, json, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
from scipy.stats import beta, norm
import sil_validation as s

ROOT=s.ROOT
OUT=ROOT/'data/independent_sil'
OUT.mkdir(exist_ok=True)

def digest(path): return hashlib.sha256(path.read_bytes()).hexdigest()

def load():
    old=s.load()
    dll=ct.CDLL(str(ROOT/'validation/sil_reporting/sil_reporting.dll'))
    for name in ('sil_control','sil_observe','sil_fit','sil_advance'):
        getattr(dll,name).argtypes=getattr(old,name).argtypes
        getattr(dll,name).restype=getattr(old,name).restype
    dll.sil_mc_run.argtypes=[ct.c_int,ct.c_int,ct.c_int,ct.c_double,s.DP,s.DP,s.DP,s.IP,s.DP]
    dll.sil_mc_run.restype=ct.c_int
    dll.sil_mc_observe.argtypes=[ct.c_int,s.DP,s.IP,s.DP,s.DP,s.DP,s.IP]
    dll.sil_mc_observe.restype=None
    dll.sil_fit_reporting.argtypes=[ct.c_int,s.DP,s.DP,s.IP,ct.c_double,s.DP,s.DP,s.DP]
    dll.sil_fit_reporting.restype=ct.c_int
    return dll

def record(dll,bits,aged,seed,N=4096,save=None):
    noise=np.random.default_rng(seed).normal(0,.002,(N+1,3))
    current=np.zeros((N+1,3));sensed=np.zeros_like(current)
    states=np.zeros((N,3),np.int32);junction=np.zeros((N+1,24))
    status=dll.sil_mc_run(N,bits,int(aged),75.,s.dp(noise),s.dp(current),s.dp(sensed),s.ip(states),s.dp(junction))
    if status:raise RuntimeError(('plant',status,bits,seed))
    X=np.zeros((N,2,50));y=np.zeros((N,2));keep=np.zeros(N,np.int32)
    dll.sil_mc_observe(N,s.dp(sensed),s.ip(states),s.dp(junction),s.dp(X),s.dp(y),s.ip(keep))
    target=np.zeros(2);cov=np.zeros((2,2));diag=np.zeros(4)
    sigma_weight=float(np.sqrt(.002**2+(64/2**bits)**2/12))
    status=dll.sil_fit_reporting(N,s.dp(X),s.dp(y),s.ip(keep),sigma_weight,s.dp(target),s.dp(cov),s.dp(diag))
    result=dict(seed=seed,aged=bool(aged),bits=bits,periods=N,status=status,
                estimate=target.tolist(),covariance=cov.tolist(),retained=int(diag[0]),
                current_sha256=hashlib.sha256(current.tobytes()).hexdigest(),
                gates_sha256=hashlib.sha256(states.tobytes()).hexdigest(),
                noise_sha256=hashlib.sha256(noise.tobytes()).hexdigest())
    if not status:
        assert np.linalg.eigvalsh(cov).min()>0 and diag[1]==36
    if save:
        np.savez_compressed(save,current=current,sensed=sensed,states=states,junction=junction,
                            X=X,y=y,keep=keep,noise=noise,target=target,covariance=cov)
    return result

def parity(dll):
    results=[]
    for bits in (16,12):
        seed=202610080001+bits
        a=record(dll,bits,True,seed,save=OUT/f'parity_adc{bits}.npz')
        z=np.load(OUT/f'parity_adc{bits}.npz')
        py=s.run(dll,4096,dict(name='independent_parity',Tcase=75,bits=bits,noise=.002),True,seed=seed)
        difference=max(float(np.max(np.abs(z[key]-py[key]))) for key in ('current','sensed','junction'))
        assert np.array_equal(z['states'],py['states']) and difference<1e-6,(bits,difference)
        est,se,rank,defect,P,_=s.p.fit(z['X'].reshape(-1,50),z['y'].ravel(),keep=z['keep'].astype(bool))
        factor=(.002**2+(64/2**bits)**2/12)/.002**2
        expected=P@P.T*factor
        np.testing.assert_allclose(a['estimate'],est,rtol=1e-9,atol=1e-10)
        np.testing.assert_allclose(a['covariance'],expected,rtol=1e-8,atol=1e-15)
        assert rank==36 and np.max(defect)<1e-7
        results.append(dict(bits=bits,max_state_difference=difference,rank=rank,
                            covariance_max_difference=float(abs(np.array(a['covariance'])-expected).max())))
    (OUT/'parity.json').write_text(json.dumps(results,indent=2)+'\n')
    print('NATIVE_PARITY',json.dumps(results),flush=True)

def rate(values):
    values=np.asarray(values,dtype=bool);n=len(values);k=int(values.sum())
    return dict(successes=k,trials=n,rate=k/n,
                exact95=[float(beta.ppf(.025,k,n-k+1)) if k else 0.,
                         float(beta.ppf(.975,k+1,n-k)) if k<n else 1.])

def difference(a,b):
    if a['status'] or b['status']:return None
    covariance=np.array(a['covariance'])+np.array(b['covariance'])
    return dict(estimate=(np.array(a['estimate'])-b['estimate']).tolist(),
                covariance=covariance.tolist(),SE=np.sqrt(np.diag(covariance)).tolist())

def summarize(rows,args):
    summary=dict(periods=args.periods,Ts_s=50e-6,calibration_pairs=args.calibration,
                 evaluation_triplets=args.evaluation,case_C=75.,hardware=False,
                 noise_A=.002,weight_rule='sigma^2 = (2 mA)^2 + LSB^2/12; MA(1) midpoint reference only',
                 pairing='Independent PCG64 streams for reference, null-test, injected-test. Common-noise comparator separately labeled.',
                 calibration_rule='Per-target order statistic of abs(null difference)/nominal independent-record SE; rank ceil(.95*(n_cal+1)); no validation/injection labels used.',
                 scope='Offline exchangeable simulated null calibration, empirical held-out alternative coverage and detection; no physical accuracy or universal GUM coverage.',scenarios={})
    for bits in (16,12):
        rr=sorted([r for r in rows if r['bits']==bits],key=lambda r:r['trial'])
        cal=[difference(r['null'],r['reference']) for r in rr if r['partition']=='calibration']
        assert all(x is not None for x in cal),'Calibration rejection must be investigated, not silently discarded'
        scores=np.array([np.abs(x['estimate'])/x['SE'] for x in cal])
        rank=int(np.ceil(.95*(len(cal)+1)));q=np.sort(scores,axis=0)[rank-1]
        ev=[r for r in rr if r['partition']=='evaluation'];positive=[];null=[];common=[];reject=0
        for r in ev:
            d=difference(r['injected'],r['reference']);n=difference(r['null'],r['reference'])
            if d is None or n is None:reject+=1;continue
            positive.append(d);null.append(n)
            if 'common_injected' in r and not r['common_injected']['status']:
                common.append(dict(common=(np.array(r['common_injected']['estimate'])-r['reference']['estimate']).tolist(),independent=d['estimate']))
        assert positive
        pe=np.array([x['estimate'] for x in positive]);ne=np.array([x['estimate'] for x in null])
        ps=np.array([x['SE'] for x in positive]);ns=np.array([x['SE'] for x in null]);truth=np.array([.05,.0002])
        scenario=dict(accepted=len(positive),rejected=reject,calibration_quantiles=q.tolist(),
                      calibration_rank=rank,mean=pe.mean(axis=0).tolist(),bias=(pe-truth).mean(axis=0).tolist(),
                      SD=pe.std(axis=0,ddof=1).tolist(),RMSE=np.sqrt(np.mean((pe-truth)**2,axis=0)).tolist(),
                      median_nominal_SE=np.median(ps,axis=0).tolist(),median_calibrated_halfwidth=np.median(q*ps,axis=0).tolist(),
                      null_mean=ne.mean(axis=0).tolist(),targets=[])
        for j in range(2):
            scenario['targets'].append(dict(name='voltage' if j==0 else 'resistance',
                nominal_coverage=rate(abs(pe[:,j]-truth[j])<=norm.ppf(.975)*ps[:,j]),
                calibrated_alternative_coverage=rate(abs(pe[:,j]-truth[j])<=q[j]*ps[:,j]),
                nominal_false_alarm=rate(abs(ne[:,j])>norm.ppf(.975)*ns[:,j]),
                calibrated_false_alarm=rate(abs(ne[:,j])>q[j]*ns[:,j]),
                nominal_detection=rate(abs(pe[:,j])>norm.ppf(.975)*ps[:,j]),
                calibrated_detection=rate(abs(pe[:,j])>q[j]*ps[:,j])))
        if common:
            common_sd=np.array([x['common'] for x in common]).std(axis=0,ddof=1)
            independent_sd=np.array([x['independent'] for x in common]).std(axis=0,ddof=1)
            scenario['common_noise_comparator']=dict(trials=len(common),common_SD=common_sd.tolist(),independent_SD=independent_sd.tolist(),SD_ratio=(common_sd/independent_sd).tolist())
        summary['scenarios'][str(bits)]=scenario
    summary['hashes']={str(p.relative_to(ROOT)).replace('\\','/'):digest(p) for p in [
        ROOT/'scripts/sil_independent_mc.py',ROOT/'scripts/sil_mc_native.c',ROOT/'scripts/sil_core.c',
        ROOT/'scripts/sil_plant.c',ROOT/'scripts/sil_basis.h',ROOT/'validation/sil_reporting/sil_reporting.dll',OUT/'trials.jsonl']}
    (ROOT/'data/independent_sil_results.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('INDEPENDENT_MC_SUMMARY',json.dumps(summary['scenarios']),flush=True)

def main(args):
    dll=load()
    if args.parity:parity(dll);return
    assert (OUT/'parity.json').is_file(),'Run native parity checks before Monte Carlo'
    path=OUT/'trials.jsonl';rows=[json.loads(line) for line in path.read_text().splitlines()] if path.exists() else []
    existing={(r['bits'],r['trial']) for r in rows}
    tasks=[(bits,trial) for bits in (16,12) for trial in range(args.calibration+args.evaluation) if (bits,trial) not in existing]
    def job(item):
        bits,trial=item;partition='calibration' if trial<args.calibration else 'evaluation'
        root_seed=202610080000+bits*1000000+trial*10
        r=dict(bits=bits,trial=trial,partition=partition)
        for name,aged,offset in [('reference',False,0),('null',False,1)]+([('injected',True,2)] if partition=='evaluation' else []):
            save=OUT/f'adc{bits}_trial{trial}_{name}.npz' if trial in (0,args.calibration,args.calibration+1) else None
            r[name]=record(dll,bits,aged,root_seed+offset,args.periods,save)
        if partition=='evaluation' and trial<args.calibration+100:
            r['common_injected']=record(dll,bits,True,root_seed,args.periods)
        assert r['reference']['noise_sha256']!=r['null']['noise_sha256']
        if 'injected' in r:assert r['reference']['noise_sha256']!=r['injected']['noise_sha256']
        return r
    start=time.perf_counter()
    with path.open('a',encoding='utf-8') as stream, ThreadPoolExecutor(max_workers=args.workers) as pool:
        for i,r in enumerate(pool.map(job,tasks)):
            stream.write(json.dumps(r)+'\n');stream.flush();rows.append(r)
            if (i+1)%20==0 or i==0:print('MC_PROGRESS',len(rows),'of',2*(args.calibration+args.evaluation),'elapsed_s',round(time.perf_counter()-start,1),flush=True)
    assert len(rows)==2*(args.calibration+args.evaluation)
    summarize(rows,args)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--parity',action='store_true');ap.add_argument('--periods',type=int,default=4096)
    ap.add_argument('--calibration',type=int,default=200);ap.add_argument('--evaluation',type=int,default=400);ap.add_argument('--workers',type=int,default=4)
    main(ap.parse_args())
