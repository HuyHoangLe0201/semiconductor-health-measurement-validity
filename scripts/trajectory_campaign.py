"""Predeclared noiseless controller ensemble: support is tested before precision.
No fixed 36-column quotient is used in analysis. Covariance is a reference,
not an empirical feedback likelihood or a statement about physical accuracy.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1'
from pathlib import Path
import ctypes as ct, itertools, json, hashlib, time, csv
import numpy as np
import sil_validation as sv
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/trajectory_campaign'; OUT.mkdir(parents=True,exist_ok=True)
DP=ct.POINTER(ct.c_double); IP=ct.POINTER(ct.c_int)
dp=lambda a:a.ctypes.data_as(DP)
ip=lambda a:a.ctypes.data_as(IP)
dll=ct.CDLL(str(ROOT/'validation/trajectory_native/trajectory.dll'))
dll.tv_set.argtypes=[ct.c_double]*6
dll.tv_set_refine.argtypes=[ct.c_int]; dll.tv_set_replay.argtypes=[ct.c_int]
dll.sil_mc_run.argtypes=[ct.c_int,ct.c_int,ct.c_int,ct.c_double,DP,DP,DP,IP,DP]
dll.sil_mc_observe.argtypes=[ct.c_int,DP,IP,DP,DP,DP,IP]
dll.sil_control.argtypes=[ct.c_double,DP,IP,IP]
dll.tv_ode.argtypes=[ct.c_double,DP,IP,ct.c_double,DP]
def run(cfg,n=4096,refine=4,replay=None):
 dll.tv_set(.002,cfg['R'],cfg['E'],cfg['A'],cfg['mod'],cfg['phase'])
 dll.tv_set_refine(refine);dll.tv_set_replay(int(replay is not None))
 current=np.zeros((n+1,3));sense=current.copy();gates=np.zeros((n,3),np.int32) if replay is None else replay.copy()
 junction=np.zeros((n+1,24));noise=np.zeros_like(current)
 status=dll.sil_mc_run(n,0,0,cfg['Tc'],dp(noise),dp(current),dp(sense),ip(gates),dp(junction))
 X=np.zeros((n,2,50));y=np.zeros((n,2));keep=np.zeros(n,np.int32)
 if status==0:dll.sil_mc_observe(n,dp(sense),ip(gates),dp(junction),dp(X),dp(y),ip(keep))
 return status,current,gates,junction,X,y,keep
scale=np.r_[np.ones(24),np.full(24,.05),1,.05]
targets=np.zeros((2,50));targets[0,[0,7]]=[1,-1];targets[1,[24,30]]=[1,-1]
ts=targets*scale
def analyze(X,keep,R,n):
 idx=np.flatnonzero(keep[:n]); rows=X[idx].reshape(-1,50)*scale
 base=dict(window=n,retained=len(idx),rows=len(rows))
 if len(rows)==0:return base|dict(rank=0,smin=None,smax=None,rank_sensitivity=[0]*3,targets=[dict(supported=False,defect=1,se=None,kappa=None) for _ in range(2)])
 _,s,V=np.linalg.svd(rows/np.sqrt(len(rows)),full_matrices=False)
 ranks=[int(np.sum(s>s[0]*tol)) for tol in (1e-8,1e-10,1e-12)]
 r=ranks[1];Vr=V[:r];coef=ts@Vr.T
 defect=np.linalg.norm(ts-coef@Vr,axis=1)/np.linalg.norm(ts,axis=1)
 kappa=s[0]*np.linalg.norm(coef/s[:r],axis=1)/np.linalg.norm(ts,axis=1)
 # Finite MA(1) covariance, restart only at a missing original period.
 c1=40+R/2;c0=-40+R/2;diag=(c1*c1+c0*c0)*.002**2;off=c1*c0*.002**2
 WX=np.empty((len(idx),2,50));prev=-2;last=0;prior=np.zeros((2,50))
 for a,k in enumerate(idx):
  sub=off/last if k==prev+1 else 0;last=np.sqrt(diag-sub*sub)
  WX[a]=(X[k]*scale-sub*prior)/last;prior=WX[a];prev=k
 _,sw,Vw=np.linalg.svd(WX.reshape(-1,50),full_matrices=False)
 rw=int(np.sum(sw>sw[0]*1e-10));cw=ts@Vw[:rw].T
 dw=np.linalg.norm(ts-cw@Vw[:rw],axis=1)/np.linalg.norm(ts,axis=1)
 se=np.linalg.norm(cw/sw[:rw],axis=1)
 result=[]
 for j in range(2):
  supported=bool(defect[j]<=1e-8 and dw[j]<=1e-8)
  result.append(dict(supported=supported,defect=float(defect[j]),whitened_defect=float(dw[j]),se=float(se[j]) if supported else None,kappa=float(kappa[j]) if supported else None))
 return base|dict(rank=r,smin=float(s[r-1]) if r else None,smax=float(s[0]),rank_sensitivity=ranks,whitened_rank=rw,targets=result)
def controller(cfg,t,sense,old):
 p=sv.p;cand=p.STATES[np.all(abs(p.STATES-old)<=1,axis=1)]
 v,r=sv.nominal(np.full(24,25));ids=p.path_ids(cand,sense,'Ttype')
 drop=np.sign(sense)*(sv.A@v[:8])[ids]+sense*(sv.A@r[:8])[ids]
 emf=cfg['E']*np.sin(2*np.pi*25*t+p.SHIFT)
 pred=sense+.025*((150*cand-drop-emf)@p.P0-cfg['R']*sense)
 nt=t+50e-6;amp=cfg['A']*(1+cfg['mod']*(5*np.sin(2*np.pi*7*nt)+3*np.sin(2*np.pi*13*nt))/16)
 ref=amp*np.sin(2*np.pi*25*nt+p.SHIFT+cfg['phase'])
 return cand[np.argmin(np.sum((pred-ref)**2,axis=1)+.005*np.sum((cand-old)**2,axis=1))]
def derivative(cfg,t,x,q):
 p=sv.p;i=x[:3];T=cfg['Tc']+x[3:].reshape(24,4).sum(axis=1);v,r=sv.nominal(T)
 active=sv.A[p.path_ids(q,i,'Ttype')]
 drop=(active*(v.reshape(3,8)+abs(i[:,None])*r.reshape(3,8))).sum(axis=1)
 pole=150*q-np.sign(i)*drop;emf=cfg['E']*np.sin(2*np.pi*25*t+p.SHIFT)
 di=(p.P0@(pole-emf)-cfg['R']*i)/.002
 power=(active*(abs(i[:,None])*v.reshape(3,8)+i[:,None]**2*r.reshape(3,8))).ravel()
 return np.r_[di,((sv.R*power[:,None]-x[3:].reshape(24,4))/sv.TAU).ravel()]
def sha(a):return hashlib.sha256(a.tobytes()).hexdigest()
def main():
 start=time.time();audit={};default=dict(R=.25,E=20,A=16,mod=1,phase=-.45,Tc=75)
 status,cur,g,T,X,y,keep=run(default,128)
 original=ct.CDLL(str(ROOT/'validation/sil_reporting/sil_reporting.dll'));original.sil_mc_run.argtypes=dll.sil_mc_run.argtypes
 cc=cur*0;ss=cur*0;gg=g*0;tt=T*0;noise=cur*0
 assert original.sil_mc_run(128,0,0,75,dp(noise),dp(cc),dp(ss),ip(gg),dp(tt))==0
 audit['default_current_max_abs']=float(np.max(abs(cur-cc)));audit['default_junction_max_abs']=float(np.max(abs(T-tt)));audit['default_gates_identical']=bool(np.array_equal(g,gg))
 assert audit['default_current_max_abs']<1e-10 and audit['default_junction_max_abs']<1e-9 and audit['default_gates_identical']
 configs=[dict(R=R,E=E,A=A,phase=phase,Tc=Tc,mod=mod) for R,E,A,phase,Tc,mod in itertools.product((.15,.4),(20,160),(2,8,20),(-.9,0,.9),(25,75),(0,1))]
 contract=dict(record_periods=4096,windows=[128,512,2048,4096],configurations=configs,Ts_s=50e-6,Vdc_V=300,L_H=.002,f_Hz=25,blanking_s=1e-6,current_exclusion_A=.5,controller_penalty=.005,acquisition='noiseless, synchronous, no quantization',thermal='cold Foster modes; oracle junction compensation',parameter_scale=scale.tolist(),rank_tolerances=[1e-8,1e-10,1e-12],support_tolerance=1e-8,reference_current_noise_A=.002,probability_interpretation='finite ensemble fractions only, not mission probabilities')
 (OUT/'contract.json').write_text(json.dumps(contract,indent=2))
 records=[];mismatch=0;count=0;derror=0;refinements=[]
 rng=np.random.default_rng(20261008)
 for k,cfg in enumerate(configs):
  status,cur,g,T,X,y,keep=run(cfg)
  entry=dict(id=k,config=cfg,status=status,current_sha256=sha(cur),gate_sha256=sha(g),junction_sha256=sha(T))
  if status==0:
   assert np.all(np.isfinite(X))
   entry['peak_current_A']=float(np.max(abs(cur)));entry['junction_range_C']=[float(T.min()),float(T.max())]
   entry['windows']=[analyze(X,keep,cfg['R'],n) for n in contract['windows']]
   for j in np.linspace(0,4095,16,dtype=int):
    q=np.zeros(3,np.int32);old=g[j-1].copy() if j else np.zeros(3,np.int32)
    dll.sil_control(j*50e-6,dp(cur[j].copy()),ip(old),ip(q));mismatch+=int(not np.array_equal(q,controller(cfg,j*50e-6,cur[j],old)));count+=1
   for _ in range(3):
    x=np.r_[rng.uniform(-20,20,3),rng.uniform(0,5,96)];q=rng.integers(-1,2,3,dtype=np.int32);t=rng.uniform(0,.2);out=np.zeros(99)
    dll.tv_ode(t,dp(x),ip(q),cfg['Tc'],dp(out));derror=max(derror,float(np.max(abs(out-derivative(cfg,t,x,q)))))
   if k in (0,23,72,143):
    s2,c2,g2,T2,*_=run(cfg,4096,8,g)
    assert s2==0 and np.array_equal(g,g2)
    refinements.append(dict(id=k,max_current_change_A=float(np.max(abs(cur-c2))),max_junction_change_C=float(np.max(abs(T-T2)))))
    np.savez_compressed(OUT/f'case_{k:03d}.npz',current=cur,gates=g,junction=T,X=X,y=y,keep=keep)
  else:entry['windows']=[]
  records.append(entry)
  if (k+1)%24==0: print(f'{k+1}/{len(configs)} complete',flush=True)
 audit.update(controller_vectors=count,controller_mismatches=mismatch,independent_derivative_max_abs=derror,fixed_gate_refinements=refinements)
 assert mismatch==0 and derror<1e-7
 assert all(w['rank']<=36 for r in records for w in r['windows'])
 summary={}
 for n in contract['windows']:
  ws=[w for r in records for w in r['windows'] if w['window']==n]
  summary[str(n)]=dict(accepted=len(ws),rank_counts={str(v):sum(w['rank']==v for w in ws) for v in sorted(set(w['rank'] for w in ws))},support_counts=[sum(w['targets'][j]['supported'] for w in ws) for j in range(2)],rank_sensitive=sum(len(set(w['rank_sensitivity']))>1 for w in ws),se_ranges=[[min([w['targets'][j]['se'] for w in ws if w['targets'][j]['supported']],default=None),max([w['targets'][j]['se'] for w in ws if w['targets'][j]['supported']],default=None)] for j in range(2)])
 result=dict(contract=contract,audit=audit,summary=summary,records=records,elapsed_s=time.time()-start)
 (OUT/'results.json').write_text(json.dumps(result,indent=2))
 (ROOT/'validation/Trajectory_campaign_QA.json').write_text(json.dumps(dict(audit=audit,summary=summary,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),status_counts={str(s):sum(r['status']==s for r in records) for s in set(r['status'] for r in records)}),indent=2))
 with (OUT/'windows.csv').open('w',newline='') as f:
  names=['id','R','E','A','phase','Tc','mod','window','retained','rank','smin','smax','target','supported','defect','whitened_defect','se','kappa'];wr=csv.DictWriter(f,fieldnames=names);wr.writeheader()
  for r in records:
   for w in r['windows']:
    for j,z in enumerate(w['targets']):wr.writerow(dict(id=r['id'],**r['config'],**{a:w[a] for a in ('window','retained','rank','smin','smax')},target=('voltage','resistance')[j],**z))
 print(json.dumps(dict(summary=summary,audit=audit,elapsed_s=time.time()-start),indent=2))
if __name__=='__main__':main()
