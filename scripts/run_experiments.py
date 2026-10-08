"""Numerical studies for the P1 MSSP manuscript.

Input designs are the audited, synthetic balanced-current operators supplied in
../data. This program does not simulate machine dynamics or a switching plant.
Run: python scripts/run_experiments.py (from any working directory).
"""
from pathlib import Path
import itertools
import json
import numpy as np
import sympy as sp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data'; FIG=ROOT/'figures'
FIG.mkdir(exist_ok=True)
RNG=np.random.default_rng(20261008)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.spines.top':False,
                     'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42,
                     'lines.linewidth':1.6,'figure.dpi':120})
COLORS=['#22577a','#c05a36','#438565']

def save(fig,name):
    fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight')
    fig.savefig(FIG/f'{name}.png',dpi=180,bbox_inches='tight')
    plt.close(fig)

def white(X,a,b,sigmaq):
    """Apply inverse lower Cholesky of tridiagonal MA(1) covariance in O(N)."""
    if a==0: return X/(b*sigmaq)
    n=len(X)//2
    d=b*b*sigmaq*sigmaq*(1+a*a)
    off=-a*b*b*sigmaq*sigmaq
    Ldiag=np.empty(n); Lsub=np.zeros(n)
    Ldiag[0]=np.sqrt(d)
    for k in range(1,n):
        Lsub[k]=off/Ldiag[k-1]
        Ldiag[k]=np.sqrt(d-Lsub[k]**2)
    shape=X.shape
    blocks=X.reshape(n,2,-1)
    result=np.empty_like(blocks)
    result[0]=blocks[0]/Ldiag[0]
    for k in range(1,n):result[k]=(blocks[k]-Lsub[k]*result[k-1])/Ldiag[k]
    return result.reshape(shape)

def estimator(X,a=.995,b=4.65,sigmaq=.0846):
    norms=np.linalg.norm(X,axis=0);norms=np.where(norms>0,norms,1)
    Xw=white(X,a,b,sigmaq)
    U,s,Vh=np.linalg.svd(Xw/norms,full_matrices=False)
    r=np.sum(s>1e-10)
    # Maps whitened observations to an arbitrary representative; targets must be estimable.
    P=(Vh[:r].T/s[:r])@U[:,:r].T/norms[:,None]
    return P,r,Xw

def ma_noise(periods,trials,a=.995,b=4.65,sigmaq=.0846):
    nu=RNG.normal(0,sigmaq,size=(periods+1,2,trials))
    return (b*(nu[1:]-a*nu[:-1])).reshape(2*periods,trials)

def rank(X):
    s=np.linalg.norm(X,axis=0);return int(np.linalg.matrix_rank(X/np.where(s>0,s,1),tol=1e-10))

def physical_mc():
    base=np.load(DATA/'Ttype_design.npz')['X']
    L=np.zeros((2,50));L[0,0]=1;L[0,7]=-1;L[1,24]=1;L[1,30]=-1
    theta=np.zeros(50);theta[0]=.05;theta[24]=.0002;theta[-2]=.15;theta[-1]=.003
    truth=L@theta
    rows=[]
    for N in [180,405,810,1620]:
        X=np.tile(base,(int(np.ceil(N/(len(base)//2))),1))[:2*N]
        P,r,Xw=estimator(X)
        # This verifies target invariance against an arbitrary representative.
        assert np.allclose(L@P@Xw,L,atol=1e-9)
        cov=L@P@P.T@L.T
        # OLS with a scalar white variance still has a calculable colored-noise covariance;
        # the deliberately naive intervals instead use sigma_w^2 (X'X)^+.
        P0,r0,X0=estimator(X,a=0,b=1,sigmaq=1)
        assert r==r0==36
        sigmaw2=4.65**2*.0846**2*(1+.995**2)
        naive_cov=sigmaw2*L@P0@P0.T@L.T
        errors=[]; naive_errors=[]
        for batch in range(10):
            noise=ma_noise(N,1000)
            # Known linear truth cancels from the estimation error exactly.
            errors.append(L@P@white(noise,.995,4.65,.0846))
            naive_errors.append(L@P0@noise)
        error=np.column_stack(errors);naive_error=np.column_stack(naive_errors)
        variance=np.var(error,axis=1,ddof=1)
        coverage=np.mean(np.abs(error)<=1.95996398454*np.sqrt(np.diag(cov))[:,None],axis=1)
        naive_coverage=np.mean(np.abs(naive_error)<=1.95996398454*np.sqrt(np.diag(naive_cov))[:,None],axis=1)
        rows.append({'periods':N,'trials':10000,'rank':int(r),
                     'truth':[float(t) for t in truth],
                     'CRLB_SE':[float(t) for t in np.sqrt(np.diag(cov))],
                     'empirical_SE':[float(t) for t in np.sqrt(variance)],
                     'variance_to_CRLB':[float(t) for t in variance/np.diag(cov)],
                     'bias':[float(t) for t in error.mean(axis=1)],
                     'coverage':[float(t) for t in coverage],
                     'naive_white_coverage':[float(t) for t in naive_coverage]})
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.9),layout='constrained')
    for j in [0,1]:
        scale=1000
        ax[j].loglog([r['periods'] for r in rows],[scale*r['CRLB_SE'][j] for r in rows],color=COLORS[j],label='Gaussian bound')
        ax[j].scatter([r['periods'] for r in rows],[scale*r['empirical_SE'][j] for r in rows],color=COLORS[j],marker='o',label='GLS: 10,000 trials')
        ax[j].set(xlabel='Observation periods',ylabel='Standard error (mV)' if j==0 else r'Standard error (m$\Omega$)',
                  title='$v_{S1a}-v_{D4a}$' if j==0 else '$r_{S1a}-r_{S4a}$')
        ax[j].grid(alpha=.2);ax[j].legend(fontsize=7)
    save(fig,'physical_precision')
    fig,ax=plt.subplots(1,2,figsize=(7.2,2.8),layout='constrained')
    for j in [0,1]:
        ax[j].plot([r['periods'] for r in rows],[100*r['coverage'][j] for r in rows],'-o',color=COLORS[0],label='GLS / MA(1) covariance')
        ax[j].plot([r['periods'] for r in rows],[100*r['naive_white_coverage'][j] for r in rows],'--s',color=COLORS[1],label='OLS / white covariance')
        ax[j].axhline(95,color='gray',lw=.8)
        ax[j].set(xlabel='Observation periods',ylabel='Coverage (%)',title='Threshold contrast' if j==0 else 'Resistance contrast')
        ax[j].legend(fontsize=7);ax[j].grid(alpha=.2)
    save(fig,'coverage')
    return rows

def nuisance_extension():
    file=np.load(DATA/'Ttype_design.npz');X=file['X'];i=file['currents_A'];n=file['switching_indicator']
    C=np.array([[1,-1,0],[1,1,-2]],float);C/=np.linalg.norm(C,axis=1)[:,None]
    Hphase=np.vstack([-C@np.diag(nn*np.sign(ii)) for nn,ii in zip(n,i)])
    assert np.allclose(Hphase.sum(axis=1),X[:,-2])
    full=np.column_stack([X[:,:48],Hphase,X[:,-1]])
    c=np.zeros(48);c[0]=1;c[7]=-1
    results={}
    for name,operator in [('shared',X),('phase_specific',full)]:
        P,r,Xw=estimator(operator)
        l=np.r_[c,np.zeros(operator.shape[1]-48)]
        assert np.allclose(l@P@Xw,l,atol=1e-9)
        variance=float(l@P@P.T@l)
        results[name]={'rank':int(r),'parameters':operator.shape[1],'threshold_contrast_SE_V':float(np.sqrt(variance))}
    return results

def shared_fixture():
    rows=[]
    for p,u,n in itertools.product(range(2),[1,2],[0,1]):
        e=np.eye(2,dtype=int)[p];rows.append(np.r_[e,u*e,n])
    X=sp.Matrix(np.array(rows).tolist());Z=X[:,:2];H=X[:,2:]
    F=Z.T*(sp.eye(8)-H*(H.T*H).inv()*H.T)*Z
    assert F==sp.Matrix([[sp.Rational(21,55),-sp.Rational(1,55)],[-sp.Rational(1,55),sp.Rational(21,55)]])
    c=sp.Matrix([1,-1]);assert (c.T*F.inv()*c)[0]==5
    # Constant switching: the absolute-level gauge survives but the contrast does not.
    fixed=np.array(rows,float);fixed[:,-1]=1
    HH=fixed[:,2:];ZZ=fixed[:,:2]
    Feff=ZZ.T@(np.eye(8)-HH@np.linalg.pinv(HH))@ZZ
    assert rank(fixed)==4
    assert np.allclose(np.array([1,-1])@np.linalg.pinv(Feff)@np.array([1,-1]),5)
    return {'exact_F':[[str(v) for v in row] for row in F.tolist()],
            'contrast_bound':5,'constant_switching_contrast_bound':5,
            'independent_path_dt_contrast_bound':5.5}

def ambiguity_figures():
    Xn=np.load(DATA/'NPC_design.npz')['X'];Xt=np.load(DATA/'Ttype_design.npz')['X']
    a=np.zeros(62);a[0]=a[8]=.05;b=np.zeros(62);b[2]=.05
    single=a.copy();single[8]=0
    t1=np.zeros(50);t1[0]=.05;t2=np.zeros(50);t2[7]=.05
    assert np.max(np.abs(Xn@(a-b)))==0
    fig,ax=plt.subplots(1,3,figsize=(7.8,2.6),layout='constrained')
    pairs=[(Xn,a,b,'NPC: multiple-device patterns'),(Xn,single,b,'NPC: single-device patterns'),(Xt,t1,t2,'T-type: outer-device patterns')]
    for panel,(X,first,second,title) in zip(ax,pairs):
        m=100;panel.plot(1000*(X@first)[:m],color=COLORS[0],label='Pattern A')
        panel.plot(1000*(X@second)[:m],'--',color=COLORS[1],label='Pattern B')
        panel.set(title=title,xlabel='Residual coordinate')
        panel.legend(fontsize=7)
    ax[0].set_ylabel('Residual (mV)')
    save(fig,'ambiguity')
    fig,axes=plt.subplots(1,3,figsize=(7.8,3.1),layout='constrained')
    names=[['S1','D1','S4','D4'],['S1','D1','S2','D2','S3','D3','S4','D4'],['S1','D1','S2','D2','S3','D3','S4','D4','Dc1','Dc2']]
    for ax,name,labels in zip(axes,['two_level','Ttype','NPC'],names):
        A=np.load(DATA/f'{name}_design.npz')['incidence'];ax.imshow(A,cmap='Blues',vmin=0,vmax=1,aspect='auto')
        ax.set_xticks(range(len(labels)),labels,rotation=70,fontsize=7)
        ax.set_yticks(range(len(A)),['P+','N+','P-','N-'] if name=='two_level' else ['P+','O+','N+','P-','O-','N-'])
        ax.set_title({'two_level':'Two-level','Ttype':'T-type','NPC':'3L-NPC'}[name])
        for row,col in zip(*np.where(A>0)):ax.text(col,row,'1',ha='center',va='center',color='white',fontsize=8)
    save(fig,'incidence')
    gamma=np.linspace(0,.995,250);d=1-8/np.pi**2;N=10000;sigma=.5
    level=np.sqrt(sigma**2/N*(2/d+gamma/(1-gamma)))*1000
    contrast=np.full_like(gamma,np.sqrt(sigma**2/N*4/d)*1000)
    fig,ax=plt.subplots(figsize=(5.6,2.8),layout='constrained')
    ax.plot(gamma,level,label='One path level',color=COLORS[1]);ax.plot(gamma,contrast,label='Path contrast',color=COLORS[0])
    ax.set(xlabel=r'Switching constancy $\gamma$',ylabel='Bound on standard error (mV)',yscale='log')
    ax.grid(alpha=.2);ax.legend(fontsize=8);save(fig,'shared_information')
    return {'NPC_multiple_output_max_V':float(np.max(np.abs(Xn@(a-b)))),
            'NPC_single_output_max_V':float(np.max(np.abs(Xn@(single-b)))),
            'Ttype_output_max_V':float(np.max(np.abs(Xt@(t1-t2))))}

def main():
    evidence={'seed':20261008,'scope':'Synthetic exogenous design; Gaussian current noise; no drive/HIL/bench validation.'}
    evidence['shared_fixture']=shared_fixture()
    evidence['ambiguity']=ambiguity_figures()
    evidence['physical_MC']=physical_mc()
    evidence['nuisance_extension']=nuisance_extension()
    ranks=[]
    for name in ['two_level','Ttype','NPC']:
        data=np.load(DATA/f'{name}_design.npz');X=data['X'];A=data['incidence'];d=A.shape[1]
        null=sp.Matrix(A.astype(int)).nullspace()
        kappa=sum(all(v[j]==0 for v in null) for j in range(d))
        ranks.append({'topology':name,'incidence_rank':int(sp.Matrix(A.astype(int)).rank()),
                      'devices':d,'kappa':int(kappa),'joint_rank':rank(X),'periods':len(X)//2})
    evidence['rank_checks']=ranks
    (DATA/'manuscript_results.json').write_text(json.dumps(evidence,indent=2),encoding='utf-8')
    # Small generated file keeps the manuscript's numerical values in sync.
    lines=['% Generated by scripts/run_experiments.py; all values are computed.']
    macros={}
    for j,tag in [(0,'Threshold'),(1,'Resistance')]:
        vals=[v['coverage'][j]*100 for v in evidence['physical_MC']]
        ratios=[v['variance_to_CRLB'][j] for v in evidence['physical_MC']]
        macros[tag+'CoverageMin']=f'{min(vals):.2f}'
        macros[tag+'CoverageMax']=f'{max(vals):.2f}'
        macros[tag+'RatioMin']=f'{min(ratios):.3f}'
        macros[tag+'RatioMax']=f'{max(ratios):.3f}'
        se=evidence['physical_MC'][-1]['CRLB_SE'][j]*1000
        macros[tag+'FinalSE']=f'{se:.3f}'
    macros['SharedContrastSE']=f"{1000*evidence['nuisance_extension']['shared']['threshold_contrast_SE_V']:.3f}"
    macros['PhaseContrastSE']=f"{1000*evidence['nuisance_extension']['phase_specific']['threshold_contrast_SE_V']:.3f}"
    for key,val in macros.items():lines.append('\\newcommand{\\'+key+'}{'+val+'}')
    (ROOT/'computed_values.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    print(json.dumps(evidence,indent=2))

if __name__=='__main__':main()
