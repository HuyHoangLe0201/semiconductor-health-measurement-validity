"""Independent mathematical and rendering audit of the monochrome figure set."""
from pathlib import Path
import json,re,hashlib,math
import numpy as np
import sympy as sp
import fitz
from scipy.stats import norm
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parents[1];F=ROOT/'figures/latex_illustrations'
fx=json.loads((F/'mathematical_fixtures.json').read_text())
A,K,h,g,zeta=[sp.Matrix(fx[k]) for k in ('A','K','h','g','zeta')]
assert A.rank()==6 and K.rank()==2 and A*K==sp.zeros(6,2)
assert A*h==-zeta and A*g==sp.ones(6,1) and K.T*K==2*sp.eye(2)
PA=A.pinv()*A
expect=sp.diag(sp.eye(4),sp.ones(2)/2,sp.ones(2)/2)
assert PA==expect and PA==sp.eye(8)-K*K.T/2
assert K.T*sp.Matrix([1,-1,0,0,0,0,0,0])==sp.zeros(2,1)
assert K.T*sp.Matrix([1,0,0,-1,0,0,0,0])==sp.zeros(2,1)
# Verify all 14 complete null directions on the paper's retained rich design.
fixture=np.load(ROOT/'data/Ttype_design.npz');X=fixture['X']
perm=[0,7,1,6,2,5,4,3]
ho=np.zeros(8);go=np.zeros(8);ko=np.zeros((8,2))
ho[perm]=np.array(h,dtype=float).ravel();go[perm]=np.array(g,dtype=float).ravel();ko[perm]=np.array(K,dtype=float)
N=np.zeros((50,14))
for block in range(6):N[8*block:8*(block+1),2*block:2*(block+1)]=ko
N[:24,12]=np.tile(ho,3);N[24:48,13]=np.tile(go,3);N[49,13]=-1
assert sp.Matrix(N.astype(int)).rank()==14
assert np.linalg.matrix_rank(X)==36
null_defect=float(abs(X@N).max());assert null_defect<1e-10,null_defect
# Exact singular example and full eigenmode relations.
Fc=sp.diag(1,sp.Rational(1,25),0);Fcp=Fc.pinv()
assert Fcp==sp.diag(1,25,0)
for kap in (0,4,36):
    mu=sp.Integer(0 if kap==0 else 1);sn=sp.Integer(1) if kap==0 else sp.Rational(1,kap)
    q=sp.ones(2,1)/2;Q=sp.eye(2)/2;d=sp.Rational(1,4)
    FF=d*Q-(mu**2*d**2/(sn+mu**2*d))*(q*q.T)
    CC=8*sp.eye(2)+kap*sp.ones(2)
    assert FF.inv()==CC
    assert CC*sp.ones(2,1)==(8+2*kap)*sp.ones(2,1)
    assert CC*sp.Matrix([1,-1])==8*sp.Matrix([1,-1])
assert sp.Matrix([[1,1,1],[1,2,2]])*sp.Matrix([0,1,-1])==sp.zeros(2,1)
# Thermal scalar fibre and linear-image vertices.
assert all(a+2*b==10 for a,b in ((0,5),(10,0),(20,-5)))
D=sp.Matrix([[1,1],[1,-1]])
vertices={tuple(D*sp.Matrix([i,j])) for i in (-1,1) for j in (-1,1)}
assert vertices=={(-2,0),(2,0),(0,-2),(0,2)}
B=sp.Matrix([[-1,1,0],[-1,0,1]])
assert B*B.T==sp.Matrix([[2,1],[1,2]])
z=norm.ppf(.975);lam=np.array(fx['power_lambda'])
pk=norm.cdf(-z-lam)+norm.sf(z-lam)
pe=norm.cdf(-z-lam/math.sqrt(2))+norm.sf(z-lam/math.sqrt(2))
np.testing.assert_allclose(pk,fx['power_known'],rtol=0,atol=1e-14)
np.testing.assert_allclose(pe,fx['power_independent_equal_records'],rtol=0,atol=1e-14)
power80=brentq(lambda t:norm.cdf(-z-t)+norm.sf(z-t)-.8,0,6)
assert abs(power80-2.801582)<1e-6 and abs(math.sqrt(2)*power80-3.962034)<1e-6
src=(F/'measurement_illustrations.tex').read_text()
curvepoints=re.findall(r'\((\d+\.\d{6}),(\d+\.\d{8})\)',src)
assert len(curvepoints)==122
for actual,want in zip(curvepoints,list(zip(lam,pk))+list(zip(lam,pe))):
    assert abs(float(actual[0])-want[0])<1e-6 and abs(float(actual[1])-want[1])<1e-7
rows=re.findall(r'^(16|12)&(Voltage|Resistance)&([\d.]+)&([\d.]+)&([\d.]+)',src,re.M)
assert len(rows)==4
data=json.loads((ROOT/'data/independent_sil_results.json').read_text())
for bits,name,q,far,detect in rows:
    j=['Voltage','Resistance'].index(name);d=data['scenarios'][bits];target=d['targets'][j]
    assert q==f"{d['calibration_quantiles'][j]:.3f}"
    assert far==f"{100*target['calibrated_false_alarm']['rate']:.2f}"
    assert detect==f"{100*target['calibrated_detection']['rate']:.2f}"
log=(F/'measurement_illustrations.log').read_text(errors='replace')
warnings=re.findall(r'(?:Overfull|Underfull|[^\n]*Warning|^!)[^\n]*',log,re.M);assert not warnings,warnings
# Exact companion scenarios; these are analytic, not new SIL outcomes.
excitation=[]
for name,u,expected in [('narrow',sp.Matrix([sp.Rational(4,5)]*2+[sp.Rational(6,5)]*2),sp.Rational(1,26)),
                       ('wide',sp.Matrix([0,0,2,2]),sp.Rational(1,2))]:
    M=sp.eye(4)-u*u.T/(u.T*u)[0]
    d0=(sp.ones(1,4)*M*sp.ones(4,1))[0]/4
    assert d0==expected and M*M==M
    var=1/(4*d0)
    excitation.append(dict(name=name,d0=str(d0),variance_over_sigma2=str(var)))
assert sp.Rational(13,2)/sp.Rational(1,2)==13
tangent=sp.Matrix([1,25,0])/sp.sqrt(26)
assert (tangent.T*Fc*tangent)[0]==1
assert sum(tangent)==sp.sqrt(26)
assert abs(20*math.exp(-math.log(4))-5)<1e-13
assert abs((1-.95**20)-.641514)<1e-6
B_reuse=sp.ones(3,1)*(-1)
B_reuse=B_reuse.row_join(sp.eye(3))
assert B_reuse*B_reuse.T==sp.Matrix([[2,1,1],[1,2,1],[1,1,2]])
alpha0=.05;zz=norm.ppf(1-alpha0/2);ss=1.7;BT=2.3
# Worst boundary biases still achieve conservative interval coverage.
for bias in np.linspace(-BT,BT,11):
    coverage=norm.cdf((zz*ss+BT-bias)/ss)-norm.cdf((-zz*ss-BT-bias)/ss)
    assert coverage>=1-alpha0-1e-12
doc=fitz.open(F/'measurement_illustrations.pdf');assert len(doc)==10
pages=[]
for i,p in enumerate(doc,1):
    assert f'Figure {i}.' in p.get_text() and 'OpenAI Codex' in p.get_text()
    assert not p.get_images()
    for drawing in p.get_drawings():
        for kind in ('color','fill'):
            color=drawing.get(kind)
            # PDF gray can acquire sub-byte RGB differences through color management.
            if color is not None:assert max(color)-min(color)<1e-3,(i,color)
    outside=[]
    for block in p.get_text('dict')['blocks']:
        for line in block.get('lines',[]):
            for s in line['spans']:
                x0,y0,x1,y1=s['bbox']
                if min(x0,y0)<-.5 or x1>p.rect.width+.5 or y1>p.rect.height+.5:outside.append(s['text'])
    assert not outside,(i,outside)
    pixels=np.frombuffer(p.get_pixmap().samples,dtype=np.uint8).reshape(-1,3)
    assert int((pixels.max(axis=1)-pixels.min(axis=1)).max())<=1
    pages.append(dict(page=i,text_within_bounds=True,vector_only=True,monochrome=True))
report=dict(edition=3,pages=10,mathematical_checks='passed',exact_incidence_kernel_and_projector=True,
 complete_design_rank=36,complete_kernel_dimension=14,complete_kernel_max_defect_V=null_defect,
 shared_information_inverse_and_eigenmodes=True,thermal_fibre_and_zonotope=True,
 estimated_baseline_covariance=True,Gaussian_power80_known=power80,
 Gaussian_power80_two_independent_records=math.sqrt(2)*power80,
 four_empirical_SIL_rows_match=True,exact_excitation_designs=excitation,excitation_information_ratio=13,
 likelihood_tangent_and_null_target=True,thermal_lag_settling_tau=math.log(4),
 bounded_bias_Gaussian_interval=True,independent_20_test_familywise_rate=1-.95**20,
 repeated_reference_covariance=True,rendering_checks=pages,LaTeX_warnings=warnings,
 PDF_SHA256=hashlib.sha256((F/'measurement_illustrations.pdf').read_bytes()).hexdigest(),
 visual_review_completed=False,hardware=False,compilation='TeX Live 2026; native compiler unavailable')
(F/'QA.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
