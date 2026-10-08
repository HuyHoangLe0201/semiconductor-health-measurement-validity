"""Classify diagnostic circuit runs and retain failed inputs/logs, never exports."""
from pathlib import Path
import json,hashlib,shutil,re

ROOT=Path(__file__).resolve().parents[1];folder=ROOT/'validation/spice';failed=ROOT/'validation/spice_failed';failed.mkdir(exist_ok=True)
records=[]
for logpath in sorted(folder.glob('*.log')):
    name=logpath.stem
    if not name.startswith(('audit_','aged75_c1_','qs_','direct_','smoothlinear_')):continue
    text=logpath.read_text(errors='replace');ended='ngspice-47 done' in text
    complete=ended and 'aborted' not in text.lower() and 'timestep too small' not in text.lower()
    time=re.search(r'Timestep too small; time = ([\d.e+-]+)',text)
    row=dict(name=name,completed=complete,scored_in_main_analysis=False,
             failure_time_s=float(time.group(1)) if time else None,
             status='completed diagnostic only' if complete else 'aborted transient' if ended else 'timeout/incomplete log',
             netlist_sha256=hashlib.sha256(logpath.with_suffix('.cir').read_bytes()).hexdigest(),log_sha256=hashlib.sha256(logpath.read_bytes()).hexdigest())
    if not complete:
        for ext in ['.cir','.log']:shutil.copy2(logpath.with_suffix(ext),failed/(name+'_failed'+ext))
    records.append(row)
result=dict(scope='Exploratory solver/model audit; no incomplete transient export is scored. Successful isolated diagnostic runs do not establish a full temperature/refinement comparison.',cases=records,
            full_window_convergence_of_original_20pF_circuit_verified=False)
(ROOT/'data/spice_extended_audit.json').write_text(json.dumps(result,indent=2))
print(json.dumps(dict(diagnostic_runs=len(records),completed=sum(r['completed'] for r in records),failed=sum(not r['completed'] for r in records)),indent=2))
