"""Independent switching-event R-L-back-EMF plant and inference diagnostics.

The plant integrates circuit dynamics; it does not synthesize y=X theta.
No parasitic capacitance, junction thermal dynamics, or hardware is represented.
"""
from pathlib import Path
import itertools, json
import numpy as np
from scipy.stats import norm
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import run_experiments as linear

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; FIG=ROOT/'figures'
C=np.array([[1,-1,0],[1,1,-2]],float)
C/=np.linalg.norm(C,axis=1)[:,None]
P0=C.T@C
STATES=np.array(list(itertools.product([-1,0,1],repeat=3)))
SHIFT=np.array([0,-2*np.pi/3,2*np.pi/3])
SEED=20261009
CFG=dict(Vdc=300.,L=.002,Rs=.25,Ts=50e-6,td=1e-6,omega=2*np.pi*25,E=20.,N=8192,sigma=.002)

def incidence(name):
    paths=([[0],[7],[1],[6]] if name=='two_level' else
           [[0],[2,5],[7],[1],[4,3],[6]] if name=='Ttype' else
           [[0,2],[8,2],[5,7],[1,3],[4,9],[4,6]])
    d=4 if name=='two_level' else 8 if name=='Ttype' else 10
    if name=='two_level': paths=[[0],[3],[1],[2]]
    A=np.zeros((len(paths),d))
    for k,path in enumerate(paths): A[k,path]=1
    return A

def path_ids(q,i,name):
    if name=='two_level': return np.where(i>=0,np.where(q==1,0,1),np.where(q==1,2,3))
    return 1-q+3*(i<0)

def emf(t,cfg): return cfg['E']*np.sin(cfg['omega']*t+SHIFT)

def reference(t):
    amp=16+5*np.sin(2*np.pi*7*t)+3*np.sin(2*np.pi*13*t)
    return amp*np.sin(CFG['omega']*t+SHIFT-.45)

def nominal_devices(d):
    v=np.where(np.arange(d)%2==0,1.1,.9)
    r=np.where(np.arange(d)%2==0,.02,.015)
    if d==10: v[-2:]=1.;r[-2:]=.015
    return v,r

def run(name='Ttype',dv=None,dr=None,cfg=None,feedback_sigma=0.,refine=2,replay=None):
    cfg=dict(CFG if cfg is None else cfg)
    A=incidence(name);d=A.shape[1];v0,r0=nominal_devices(d)
    v=np.tile(v0,(3,1));r=np.tile(r0,(3,1))
    if dv is not None: v+=dv
    if dr is not None: r+=dr
    path_v=v@A.T;path_r=r@A.T;phase_indices=np.arange(3)
    N=cfg['N'];dt=cfg['Ts'];rng=np.random.default_rng(SEED)
    feedback_noise=rng.normal(0,feedback_sigma,(N+1,3))@P0
    states=STATES if name!='two_level' else np.array(list(itertools.product([-1,1],repeat=3)))
    current=np.zeros((N+1,3));chosen=np.zeros((N,3),int);means=np.zeros((N,3))
    qold=np.zeros(3,int) if name!='two_level' else -np.ones(3,int)
    for k in range(N):
        t=k*dt;i=current[k]
        if replay is None:
            sense=i+feedback_noise[k]
            allowed=np.all(np.abs(states-qold)<=1,axis=1) if name!='two_level' else np.ones(len(states),bool)
            cand=states[allowed]
            ip=np.broadcast_to(sense,(len(cand),3));ids=path_ids(cand,ip,name)
            drop=np.sign(ip)*(A@v0)[ids]+ip*(A@r0)[ids]
            u=(cfg['Vdc']/2)*cand-drop
            pred=ip+dt/CFG['L']*((u-emf(t,CFG))@P0-CFG['Rs']*ip)
            cost=np.sum((pred-reference(t+dt))**2,axis=1)+.005*np.sum((cand-qold)**2,axis=1)
            qnew=cand[np.argmin(cost)]
        else: qnew=replay[k]
        chosen[k]=qnew
        # During all-off blanking, current freewheels toward the lower/upper level.
        qfree=np.where(i>=0,np.minimum(qold,qnew),np.maximum(qold,qnew))
        qfree=np.where(qold==qnew,qnew,qfree)
        acc=np.zeros(3);j=i.copy()
        def derivative(tt,jj,qq):
            pp=path_ids(qq,jj,name)
            pv=path_v[phase_indices,pp]
            pr=path_r[phase_indices,pp]
            pole=cfg['Vdc']/2*qq-np.sign(jj)*pv-jj*pr
            return (P0@(pole-emf(tt,cfg))-cfg['Rs']*jj)/cfg['L']
        def advance(tt,jj,qq,h,depth=0):
            k1=derivative(tt,jj,qq);s2=jj+h*k1/2;k2=derivative(tt+h/2,s2,qq)
            s3=jj+h*k2/2;k3=derivative(tt+h/2,s3,qq);s4=jj+h*k3;k4=derivative(tt+h,s4,qq)
            end=jj+h*(k1+2*k2+2*k3+k4)/6
            signs=np.sign(np.stack([jj,s2,s3,s4,end]))
            if depth<10 and np.any(signs.min(axis=0)!=signs.max(axis=0)):
                first,acc1=advance(tt,jj,qq,h/2,depth+1)
                second,acc2=advance(tt+h/2,first,qq,h/2,depth+1)
                return second,acc1+acc2
            return end,h*(jj+2*s2+2*s3+s4)/6
        offset=0.
        for qq,duration in [(qfree,cfg['td']),(qnew,dt-cfg['td'])]:
            count=1 if duration==cfg['td'] else refine
            h=duration/count
            for step in range(count):
                j,increment=advance(t+offset,j,qq,h)
                acc+=increment
                offset+=h
        current[k+1]=P0@j;means[k]=acc/dt;qold=qnew
    return dict(current=current,states=chosen,mean_current=means,cfg=cfg,name=name,feedback_noise=feedback_noise)

def observations(record,sigma=0.,seed=SEED+1,extended=False,Lref=None,use_feedback_record=False):
    name=record['name'];cfg=record['cfg'];A=incidence(name);d=A.shape[1]
    rng=np.random.default_rng(seed)
    sensed=record['current']+(record['feedback_noise'] if use_feedback_record else rng.normal(0,sigma,record['current'].shape)@P0)
    mid=(sensed[:-1]+sensed[1:])/2
    q=record['states'];old=np.vstack([np.zeros(3,int),q[:-1]])
    if name=='two_level': old[0]=-1
    # Estimator occupation uses declared gate blanking, not the plant integration.
    free=np.where(mid>=0,np.minimum(old,q),np.maximum(old,q))
    frac=CFG['td']/CFG['Ts'];pv,pr=nominal_devices(d)
    times=(np.arange(len(q))+.5)*CFG['Ts']
    e=CFG['E']*np.sin(CFG['omega']*times[:,None]+SHIFT)
    ids=path_ids(q,mid,name);fids=path_ids(free,mid,name)
    Bleg=(1-frac)*A[ids]+frac*A[fids]
    B=np.zeros((len(q),3,3*d))
    for phase in range(3): B[:,phase,phase*d:(phase+1)*d]=Bleg[:,phase]
    z=np.sign(mid);n=np.abs(q-free);delta=sensed[1:]-sensed[:-1]
    V=-np.einsum('ab,kb,kbd->kad',C,z,B)
    R=-np.einsum('ab,kb,kbd->kad',C,mid,B)
    H=np.stack([-(n*z)@C.T,-mid@C.T],axis=2)
    if extended is True:H=np.concatenate([H,(-delta@C.T/CFG['Ts'])[:,:,None],(-e@C.T/CFG['E'])[:,:,None]],axis=2)
    elif extended=='emf':H=np.concatenate([H,(-e@C.T/CFG['E'])[:,:,None]],axis=2)
    X=np.concatenate([V,R,H],axis=2).reshape(2*len(q),-1)
    nominal=CFG['Vdc']/2*q-z*(Bleg@pv)-mid*(Bleg@pr)-CFG['Rs']*mid-e
    res=((CFG['L'] if Lref is None else Lref)/CFG['Ts']*delta-nominal)@C.T
    valid=np.all((np.abs(sensed[:-1])>.5)&(np.abs(sensed[1:])>.5)&(sensed[:-1]*sensed[1:]>0),axis=1)
    return X,res.ravel(),valid

def fit(X,y,d=8,sigma=.002,keep=None,Lref=None):
    # Midpoint residual has two sensor coefficients b +/- R0/2.
    c1=(CFG['L'] if Lref is None else Lref)/CFG['Ts']+CFG['Rs']/2
    c0=-(CFG['L'] if Lref is None else Lref)/CFG['Ts']+CFG['Rs']/2
    idx=np.arange(len(y)//2) if keep is None else np.flatnonzero(keep)
    X=X.reshape(-1,2,X.shape[1])[idx].reshape(-1,X.shape[1]);y=y.reshape(-1,2)[idx].ravel()
    def whiten(values):
        shape=values.shape;blocks=values.reshape(len(idx),2,-1);out=np.empty_like(blocks)
        diagonal=(c1*c1+c0*c0)*sigma*sigma;off=c1*c0*sigma*sigma
        last=np.sqrt(diagonal);out[0]=blocks[0]/last
        for k in range(1,len(idx)):
            sub=off/last if idx[k]==idx[k-1]+1 else 0.
            last=np.sqrt(diagonal-sub*sub);out[k]=(blocks[k]-sub*out[k-1])/last
        return out.reshape(shape)
    Xw=whiten(X);norms=np.linalg.norm(X,axis=0);norms=np.where(norms>0,norms,1)
    U,s,Vh=np.linalg.svd(Xw/norms,full_matrices=False);r=int(np.sum(s>1e-10))
    P=(Vh[:r].T/s[:r])@U[:,:r].T/norms[:,None]
    l=np.zeros((2,X.shape[1]));l[0,0]=1;l[0,7]=-1;l[1,3*d]=1;l[1,3*d+6]=-1
    defect=np.max(np.abs(l@P@Xw-l),axis=1)
    theta=P@whiten(y[:,None]);estimate=l@theta
    se=np.sqrt(np.diag(l@P@P.T@l.T))
    details={'retained_periods':len(idx),'minimum_scaled_singular_value':float(s[r-1]),'representative_nuisance':theta[6*d:,0].tolist()}
    return estimate[:,0],se,int(r),defect,l@P,details

def save(fig,name):
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight');fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight');plt.close(fig)

def make_detection():
    base=json.loads((DATA/'manuscript_results.json').read_text())
    recs=base['physical_MC'];z=norm.ppf(.975)
    out=[];fig,ax=plt.subplots(1,2,figsize=(7.2,2.9),layout='constrained')
    for j in range(2):
        effects=np.linspace(0,.15 if j==0 else .001,240)
        for rec in recs:
            se=rec['CRLB_SE'][j];power=norm.sf(z-effects/se)+norm.cdf(-z-effects/se)
            ax[j].plot(1000*effects,100*power,label=f"N={rec['periods']}")
        ax[j].axhline(80,lw=.7,color='gray',ls=':');ax[j].set(xlabel='Contrast change (mV)' if j==0 else r'Contrast change (m$\Omega$)',ylabel='Detection power (%)',ylim=(0,102),title='Ideal Gaussian threshold test' if j==0 else 'Ideal Gaussian resistance test');ax[j].legend(fontsize=7);ax[j].grid(alpha=.2)
        delta=recs[-1]['truth'][j];se=recs[-1]['CRLB_SE'][j]
        out.append(dict(target=j,power=float(norm.sf(z-delta/se)+norm.cdf(-z-delta/se)),halfwidth95=float(z*se),approx_MDE80=float((z+norm.ppf(.8))*se)))
    save(fig,'detection_power');return out

def separating_fixture():
    A=incidence('Ttype');d=8;rows=[];graph=np.zeros((18,18),bool)
    signs=np.array([1,1,1,-1,-1,-1]);samples=0
    for ids in itertools.product(range(6),repeat=3):
        z=signs[list(ids)]
        if abs(z.sum())==3:continue
        minority=np.where(z!=np.sign(z.sum()))[0][0];major=[j for j in range(3) if j!=minority]
        for m1,m2 in [(40,60),(50,60),(40,70)]:
            mags=np.zeros(3);mags[major]=[m1,m2];mags[minority]=m1+m2;i=mags*z
            B=np.zeros((3,24))
            for j,p in enumerate(ids):B[j,j*d:(j+1)*d]=A[p]
            for n in [np.zeros(3),np.array([1,0,0])]:
                rows.append(np.column_stack([-C@np.diag(z)@B,-C@np.diag(i)@B,-C@(n*z),-C@i]));samples+=1
        vertices=[j*6+ids[j] for j in range(3)]
        for x in vertices:
            for y in vertices:graph[x,y]=True
    seen={0}
    while True:
        expanded=seen|set(np.where(graph[list(seen)].any(axis=0))[0])
        if expanded==seen:break
        seen=expanded
    X=np.vstack(rows);r=linear.rank(X)
    assert len(seen)==18 and r==36
    np.savez_compressed(DATA/'separating_fixture.npz',X=X)
    return dict(periods=samples,connected_vertices=len(seen),rank=r,affine_rank_per_triple=3,scope='Abstract finite separating experiment, not a plant trajectory.')

def main():
    dv=np.zeros((3,8));dr=dv.copy();dv[0,0]=.05;dr[0,0]=.0002
    record=run(dv=dv,dr=dr)
    X,y,keep=observations(record);estimate,se,r,defect,T,details=fit(X,y,keep=keep)
    assert np.max(defect)<1e-7
    print('Matched plant:',estimate,se,r,flush=True)
    output={'seed':SEED,'plant_parameters':CFG,'plant_scope':'Switching-event R-L-back-EMF simulation with imposed rotor speed, ideal devices and blanking; no parasitic, mechanical or thermal states.','fixture':separating_fixture(),'detection':make_detection(),'conditions':[]}
    runs=[('matched',record,False,None),('noisy_feedback',run(dv=dv,dr=dr,feedback_sigma=CFG['sigma']),False,None)]
    wrong=dict(CFG,L=CFG['L']*1.01,E=CFG['E']*1.01)
    mismatched=run(dv=dv,dr=dr,cfg=wrong)
    runs.extend([('L_and_EMF_mismatch',mismatched,False,None),('extended_nuisance',mismatched,True,None),('independent_L_reference',mismatched,'emf',wrong['L'])])
    for name,rec,extended,Lref in runs:
        XX,yy,valid=observations(rec,sigma=CFG['sigma'] if name=='noisy_feedback' else 0.,extended=extended,Lref=Lref,use_feedback_record=name=='noisy_feedback')
        est,ss,rr,dd,_,info=fit(XX,yy,keep=valid,Lref=Lref)
        output['conditions'].append(dict(name=name,estimate=est.tolist(),nominal_SE=ss.tolist(),rank=rr,estimability_defect=dd.tolist(),bias=(est-np.array([.05,.0002])).tolist(),**info))
        print(name,est,ss,rr,flush=True)
    # Monte Carlo regenerates complete sensor records and refits noisy regressors.
    # Gates are held to one recorded trajectory; this is not feedback Monte Carlo.
    trial_results=[];trial_errors=[]
    for trial in range(100):
        XX,yy,valid=observations(record,sigma=CFG['sigma'],seed=SEED+100+trial)
        est,ss,rr,dd,_,_=fit(XX,yy,keep=valid)
        assert rr==36 and np.max(dd)<1e-7
        trial_results.append(est)
        trial_errors.append(ss)
        if (trial+1)%25==0: print('Measurement trials:',trial+1,flush=True)
    trial_results=np.array(trial_results);trial_errors=np.array(trial_errors)
    output['measurement_MC']={'trials':len(trial_results),'scope':'Repeated sensor perturbations conditional on recorded gates; no controller rerun. Intervals use each refitted measured-design nominal SE.','mean':trial_results.mean(axis=0).tolist(),'SD':trial_results.std(axis=0,ddof=1).tolist(),'RMSE':np.sqrt(np.mean((trial_results-[.05,.0002])**2,axis=0)).tolist(),'nominal_SE':se.tolist(),'mean_refitted_SE':trial_errors.mean(axis=0).tolist(),'nominal_coverage':np.mean(np.abs(trial_results-[.05,.0002])<=1.95996398454*trial_errors,axis=0).tolist()}
    # Numerical integration convergence uses the exact same gate sequence.
    fine=run(dv=dv,dr=dr,refine=4,replay=record['states'])
    output['integration_max_difference_A']=float(np.max(np.abs(record['current']-fine['current'])))
    npc_dv=np.zeros((3,10));npc_dv[0,[0,8]]=.05
    npc_a=run('NPC',dv=npc_dv)
    npc_dv_b=np.zeros_like(npc_dv);npc_dv_b[0,2]=.05
    npc_b=run('NPC',dv=npc_dv_b,replay=npc_a['states'])
    npc_single=np.zeros_like(npc_dv);npc_single[0,0]=.05
    npc_c=run('NPC',dv=npc_single,replay=npc_a['states'])
    output['NPC_gate_replay']={'ambiguous_max_current_difference_A':float(np.max(np.abs(npc_a['current']-npc_b['current']))),'single_device_max_current_difference_A':float(np.max(np.abs(npc_c['current']-npc_b['current'])))}
    np.savez_compressed(DATA/'plant_trace.npz',current=record['current'],states=record['states'],mean_current=record['mean_current'],X=X,y=y,measurement_MC=trial_results,measurement_MC_SE=trial_errors)
    fig,ax=plt.subplots(2,1,figsize=(7.2,3.8),layout='constrained');t=np.arange(CFG['N']+1)*CFG['Ts']
    ax[0].plot(t,record['current'][:,0],label='Phase-a plant current');ax[0].plot(t,np.array([reference(tt)[0] for tt in t]),'--',lw=.9,label='Reference');ax[0].set(ylabel='Current (A)');ax[0].legend(fontsize=7)
    ax[1].step(t[:200],record['states'][:200,0],where='post');ax[1].set(xlabel='Time (s)',ylabel='Applied level',yticks=[-1,0,1]);save(fig,'plant_trace')
    fig,axes=plt.subplots(1,2,figsize=(7.2,3),layout='constrained')
    shown=[v for v in output['conditions'] if v['name']!='extended_nuisance']
    labels=[v['name'].replace('_','\n') for v in shown]
    for j,ax in enumerate(axes):
        ax.errorbar(np.arange(len(labels)),[1000*v['estimate'][j] for v in shown],yerr=[1.96*1000*v['nominal_SE'][j] for v in shown],fmt='o',capsize=3)
        ax.axhline(50 if j==0 else .2,color='gray',ls='--',label='Injected contrast');ax.set_xticks(range(len(labels)),labels,fontsize=7);ax.set(ylabel='Threshold contrast (mV)' if j==0 else r'Resistance contrast (m$\Omega$)');ax.legend(fontsize=7)
    save(fig,'plant_bias')
    (DATA/'plant_results.json').write_text(json.dumps(output,indent=2),encoding='utf-8')
    print(json.dumps(output,indent=2))

if __name__=='__main__':main()
