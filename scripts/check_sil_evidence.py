"""Independent artifact/provenance audit, read-only except the audit report."""
from pathlib import Path
import hashlib,json
import numpy as np
import fitz

ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'data/sil_results.json').read_text());assert r['hardware_validation_performed'] is False
for name,digest in r['hashes'].items():assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
source=ROOT/'external/datasheets/Infineon_IKx40N65H5.pdf';provenance=json.loads(source.with_suffix('.json').read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest()==provenance['sha256']
doc=fitz.open(source);thermal=doc[11].get_text();static=doc[4].get_text()
for val in ['0.08245484','0.144197','0.2151774','0.1581708','0.6701584','0.775759','0.3540826','0.01235548','0.08020881','0.04680901']:assert val in thermal,val
for val in ['1.65','1.85','1.95','1.45','1.40']:assert val in static,val
traces=[]
for row in r['conditions'][:8]:
    for aged in ['healthy','aged']:
        path=ROOT/'data'/f'sil_{row["name"]}_{aged}.npz';z=np.load(path);N=r['periods']
        assert z['current'].shape==(N+1,3) and z['junction'].shape==(N+1,24) and z['states'].shape==(N,3)
        for key in z.files:assert np.isfinite(z[key]).all()
        assert np.isin(z['states'],[-1,0,1]).all()
        assert abs(z['current'].sum(axis=1)).max()<1e-8
        assert z['junction'].min()>=row['Tcase']-1e-7 and z['junction'].max()<175
        if row.get('bits',0):
            scaled=z['sensed']/(64/2**row['bits']);np.testing.assert_allclose(scaled,np.rint(scaled),atol=1e-12,rtol=0)
        traces.append(dict(name=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),peak_current_A=float(abs(z['current']).max()),junction_max_C=float(z['junction'].max())))
assert r['checks']['controller_mismatches']==0 and r['checks']['nonseparating_status']<0
assert r['checks']['operator_retention_match'] and r['checks']['operator_X_difference']<1e-12 and r['checks']['operator_y_difference']<1e-11
assert r['feedback_MC']['distinct_aged_gate_records']==5
assert r['refinement']['periods']==r['periods']
coarse=np.load(ROOT/'data/sil_clean75_aged.npz');fine=np.load(ROOT/'data/sil_clean75_aged_refine8.npz')
np.testing.assert_array_equal(coarse['states'],fine['states'])
assert fine['current'].shape==(r['periods']+1,3)
assert abs(float(abs(fine['current']-coarse['current']).max())-r['refinement']['max_current_difference_A'])<1e-12
assert abs(float(abs(fine['junction']-coarse['junction']).max())-r['refinement']['max_temperature_difference_K'])<1e-12
report=dict(passed=True,source_parameters_verified_from_local_manufacturer_PDF=True,source_hash=provenance['sha256'],
 controller_operator_GLS_checks=r['checks'],trace_checks=traces,hardware_validation_performed=False)
(ROOT/'validation/sil_evidence_checks.json').write_text(json.dumps(report,indent=2))
print(json.dumps(dict(passed=True,traces=len(traces),source_hash=provenance['sha256']),indent=2))
