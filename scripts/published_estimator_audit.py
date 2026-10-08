"""Local averaged-model audit of Ou et al., APEC 2025, equations (1)--(10).
This is not a reproduced PI/PWM controller, HIL trace, or published-error benchmark.
The nominal duty is frozen for a first-order matched-operating-point derivative.
"""
from pathlib import Path
import json,itertools,numpy as np
ROOT=Path(__file__).resolve().parents[1]
out=ROOT/'data/published_estimator_audit';out.mkdir(parents=True,exist_ok=True)
rows=[];P=np.eye(3)-np.ones((3,3))/3
for amp,phase in itertools.product((10,20,40),(-np.arccos(.8),0,np.arccos(.8))):
 theta=2*np.pi*np.arange(400)/400;ang=theta[:,None]+np.array([0,-2*np.pi/3,2*np.pi/3])
 i=amp*np.sin(ang+phase);di=amp*2*np.pi*50*np.cos(ang+phase)
 vstar=325*np.sin(ang)+.1*i+.006*di;duty=.5+vstar/800
 assert duty.min()>0 and duty.max()<1
 mask=i[:,0]>=0;g=i[:,0].clip(0)*duty[:,0];den=g.sum();weights=1.5*mask/den
 rvec=np.zeros_like(i);rvec[:,0]=g
 vvec=np.zeros_like(i);vvec[:,0]=mask*duty[:,0]
 # Disturbances to required compensating phase voltage, then neutral projection.
 H=np.array([weights@((a@P)[:,0]) for a in (rvec,vvec,i,np.sign(i))])
 assert abs(H[0]-1)<1e-12
 # H acts on [dR_S1 (ohm), dV_S1 (volt), dR_filter (ohm),
 # common all-conduction-element dV (volt)]. One averaged scalar per point.
 Sc=np.array([.001,.02,.001,.01]);l=np.array([.001,0,0,0]);hs=H*Sc
 defect=float(np.linalg.norm(l-(l@hs)/(hs@hs)*hs)/np.linalg.norm(l))
 rows.append(dict(amplitude_A=amp,current_phase_rad=float(phase),power_factor=.8 if phase else 1.,duty_min=float(duty.min()),duty_max=float(duty.max()),local_operator=H.tolist(),scaled_target_defect=defect,
                  resistance_only_estimate_Ohm=.001,
                  threshold_only_apparent_Ohm=float(H[1]*.02),
                  filter_only_apparent_Ohm=float(H[2]*.001),
                  common_drop_only_apparent_Ohm=float(H[3]*.01)))
result=dict(source=dict(doi='10.1109/APEC48143.2025.10977401',equations='(1)--(10)',accepted_manuscript='https://vbn.aau.dk/ws/portalfiles/portal/779663726/2024360959.pdf'),
 scope='local noise-free averaged observation model at matched operating points; frozen nominal duty; not a controller/HIL reproduction',
 assumptions=dict(Vdc_V=800,grid_peak_V=325,Lfilter_H=.006,Rfilter_Ohm=.1,f_Hz=50,samples_per_cycle=400),
 targets_and_scales=['1 mOhm S1 resistance','20 mV S1 threshold','1 mOhm filter resistance','10 mV all-element drop'],
 baseline_cancellation='stationary nuisance cancels; changing nuisance between records does not',
 rows=rows)
(out/'results.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(max_restricted_error=max(abs(r['resistance_only_estimate_Ohm']-.001) for r in rows),
 threshold_bias_mOhm=[1000*min(r['threshold_only_apparent_Ohm'] for r in rows),1000*max(r['threshold_only_apparent_Ohm'] for r in rows)],
 filter_bias_mOhm=[1000*min(r['filter_only_apparent_Ohm'] for r in rows),1000*max(r['filter_only_apparent_Ohm'] for r in rows)],
 common_drop_bias_mOhm=[1000*min(r['common_drop_only_apparent_Ohm'] for r in rows),1000*max(r['common_drop_only_apparent_Ohm'] for r in rows)],
 defect_range=[min(r['scaled_target_defect'] for r in rows),max(r['scaled_target_defect'] for r in rows)]),indent=2))
