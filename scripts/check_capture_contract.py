"""Meaningful rejection tests for the exported capture-analysis prototype."""
from pathlib import Path
import csv,hashlib,json,shutil,tempfile
import numpy as np
import audit_capture as audit

ROOT=Path(__file__).resolve().parents[1]

def main():
    original=ROOT/'validation/offline_hot_aged.csv'
    baseline=ROOT/'validation/offline_cold_healthy.csv'
    expected=json.loads((ROOT/'data/thermal_results.json').read_text())['conditions'][3]['estimate']
    result=audit.analyze(original,baseline)
    assert result['hardware_validation_performed'] is False
    assert np.allclose(result['estimate'],expected,rtol=0,atol=1e-10)
    checks={'offline_round_trip_matches_computed_oracle':True}
    with tempfile.TemporaryDirectory(prefix='p1_capture_checks_') as temp:
        dest=Path(temp)/'record.csv';meta_path=dest.with_suffix('.json')
        def reset():
            shutil.copy2(original,dest)
            return json.loads(original.with_suffix('.json').read_text())
        def rejected(label,metadata):
            meta_path.write_text(json.dumps(metadata),encoding='utf-8')
            try:audit.read_capture(dest)
            except ValueError:checks[label]=True
            else:raise AssertionError('Invalid capture accepted: '+label)
        meta=reset();meta['evidence_kind']='controller_hil'
        rejected('offline_data_cannot_be_relabeled_as_HIL',meta)
        meta=reset();meta['csv_sha256']='0'*64
        rejected('digest_mismatch_rejected',meta)
        meta=reset();meta['temperature_basis']='case_proxy'
        rejected('case_proxy_cannot_be_promoted_to_junction',meta)
        meta=reset();meta['independent_machine_scale_reference']=False
        rejected('missing_machine_scale_reference_rejected',meta)
        meta=reset();meta['nominal_current_noise_SD_A']=.02
        rejected('unadapted_sensor_covariance_rejected',meta)
        meta=reset();meta.update(evidence_kind='controller_hil',source_method='hardware_capture',real_time_execution=True,
                                hardware_id='synthetic_rejection_test_no_hardware',capture_session_id='synthetic',independent_reference_id='synthetic')
        with dest.open(newline='',encoding='utf-8') as handle:rows=list(csv.reader(handle))
        rows[0].append('controller_latency_s')
        for row in rows[1:]:row.append('.00006')
        with dest.open('w',newline='',encoding='utf-8') as handle:csv.writer(handle).writerows(rows)
        meta['csv_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest()
        rejected('synthetic_missed_deadline_rejected',meta)
        meta=reset()
        with dest.open(newline='',encoding='utf-8') as handle:rows=list(csv.reader(handle))
        rows.pop(42)
        with dest.open('w',newline='',encoding='utf-8') as handle:csv.writer(handle).writerows(rows)
        meta['csv_sha256']=hashlib.sha256(dest.read_bytes()).hexdigest()
        rejected('missing_frame_with_valid_digest_rejected',meta)
    (ROOT/'validation/capture_contract_checks.json').write_text(json.dumps(checks,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(checks,indent=2))

if __name__=='__main__':main()
