"""Audit disjoint streams, frozen calibration, and every published rate."""
from pathlib import Path
import json, hashlib, math
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1]
raw=ROOT/'data/independent_sil/trials.jsonl'
rows=[json.loads(x) for x in raw.read_text().splitlines()]
report=json.loads((ROOT/'data/independent_sil_results.json').read_text())
assert len(rows)==1200 and len({(r['bits'],r['trial']) for r in rows})==1200
parity=json.loads((ROOT/'data/independent_sil/parity.json').read_text())
assert len(parity)==2 and all(p['rank']==36 and p['max_state_difference']<1e-6 and p['covariance_max_difference']<1e-14 for p in parity)
seeds=set();hashes=set();records=0;checks=[]
for bits in (16,12):
    rr=sorted([r for r in rows if r['bits']==bits],key=lambda r:r['trial'])
    assert [r['trial'] for r in rr]==list(range(600))
    cal=[];pos=[];null=[]
    for row in rr:
        assert row['partition']==('calibration' if row['trial']<200 else 'evaluation')
        assert ('injected' in row)==(row['trial']>=200)
        for name in ('reference','null','injected'):
            if name not in row:continue
            rec=row[name];records+=1
            assert rec['status']==0 and rec['retained']>=18 and rec['periods']==4096
            assert rec['seed'] not in seeds and rec['noise_sha256'] not in hashes
            seeds.add(rec['seed']);hashes.add(rec['noise_sha256'])
            covariance=np.array(rec['covariance'])
            assert np.allclose(covariance,covariance.T) and np.linalg.eigvalsh(covariance).min()>0
        if 'common_injected' in row:
            assert row['common_injected']['seed']==row['reference']['seed']
            assert row['common_injected']['noise_sha256']==row['reference']['noise_sha256']
            assert row['common_injected']['aged']
        def delta(name):
            d=np.array(row[name]['estimate'])-row['reference']['estimate']
            cov=np.array(row[name]['covariance'])+np.array(row['reference']['covariance'])
            return d,np.sqrt(np.diag(cov))
        n,ns=delta('null')
        if row['partition']=='calibration':cal.append(abs(n)/ns)
        else:
            p,ps=delta('injected');null.append((n,ns));pos.append((p,ps))
    rank=math.ceil(.95*201);q=np.sort(cal,axis=0)[rank-1]
    np.testing.assert_array_equal(q,report['scenarios'][str(bits)]['calibration_quantiles'])
    for j in range(2):
        pe=np.array([a[0][j] for a in pos]);ps=np.array([a[1][j] for a in pos])
        ne=np.array([a[0][j] for a in null]);ns=np.array([a[1][j] for a in null])
        truth=[.05,.0002][j]
        counts=dict(nominal_coverage=int((abs(pe-truth)<=norm.ppf(.975)*ps).sum()),
                    calibrated_alternative_coverage=int((abs(pe-truth)<=q[j]*ps).sum()),
                    nominal_false_alarm=int((abs(ne)>norm.ppf(.975)*ns).sum()),
                    calibrated_false_alarm=int((abs(ne)>q[j]*ns).sum()),
                    nominal_detection=int((abs(pe)>norm.ppf(.975)*ps).sum()),
                    calibrated_detection=int((abs(pe)>q[j]*ps).sum()))
        saved=report['scenarios'][str(bits)]['targets'][j]
        assert all(saved[k]['successes']==v and saved[k]['trials']==400 for k,v in counts.items())
        checks.append(dict(bits=bits,target=j,counts=counts))
for rel,h in report['hashes'].items():
    assert hashlib.sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
result=dict(groups=1200,primary_records=records,common_noise_comparators=200,
            distinct_primary_noise_streams=len(hashes),null_calibration_disjoint=True,
            rates_recomputed_from_trial_records=True,native_covariance_and_trajectory_parity_verified=True,
            checks=checks,hardware=False)
(ROOT/'validation/independent_sil_audit.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
