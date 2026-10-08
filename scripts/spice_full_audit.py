"""Full-window solver audit. Failed exports are never scored."""
from pathlib import Path
import json, subprocess, concurrent.futures, shutil, argparse, re
import numpy as np
import spice_validation as s

ROOT=s.ROOT
def candidate(tag, changes):
    states=np.load(ROOT/'data/thermal_cold_healthy.npz')['states'][:1024]
    path=s.write_netlist(states,tag,75,True,125e-9)
    net=path.read_text()
    for old,new in changes:
        net=re.sub(old[3:],new,net,flags=re.M) if old.startswith('RE:') else net.replace(old,new)
    path.write_text(net,encoding='ascii')
    try:r=subprocess.run([str(s.EXE),'-b','-o',tag+'.log',path.name],cwd=s.OUT,capture_output=True,text=True,timeout=1200)
    except subprocess.TimeoutExpired:
        r=subprocess.CompletedProcess([],124)
    log=path.with_suffix('.log').read_text(errors='replace')
    valid=r.returncode==0 and 'aborted' not in log.lower() and 'timestep too small' not in log.lower()
    result=dict(name=tag,completed=valid,changes=changes,returncode=r.returncode)
    if valid:
        a=np.loadtxt(s.OUT/(tag+'_current.txt'),skiprows=1)
        assert a.shape==(1025,4)
        result['neutral_sum_A']=float(abs(a[:,1:].sum(axis=1)).max())
    else:
        fail=ROOT/'validation/spice_failed';fail.mkdir(exist_ok=True)
        for ext in ['.cir','.log']:shutil.copy2(path.with_suffix(ext),fail/(tag+'_failed'+ext))
        result['failure_tail']=log[-500:]
    print(json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--second',action='store_true');parser.add_argument('--third',action='store_true');parser.add_argument('--fourth',action='store_true');parser.add_argument('--shunt',action='store_true');args=parser.parse_args()
    candidates=[
        ('audit_trap125',[('method=gear maxord=2','method=trap maxord=2')]),
        ('audit_be125',[('method=gear maxord=2','method=gear maxord=1')]),
        ('audit_esr125',[('Vpos pos 0 150','Vpos vp 0 150\nRbusp vp pos 0.01'),('Vneg neg 0 -150','Vneg vn 0 -150\nRbusn vn neg 0.01')])]
    if args.second:
        candidates=[('audit_nocap125',[('Cjo=20p','Cjo=0')]),
          ('audit_switch125',[('Cjo=20p','Cjo=0'),('RE:^Bsw([abc])([1-4]) (\\S+) (\\S+) I=.*$',r'Ssw\1\2 \3 \4 g\1\2 0 gateSW')])]
    if args.third:
        candidates=[('audit_smp125',[('set klu','unset klu')]),
                    ('audit_abstol125',[('abstol=1e-9','abstol=1e-7')])]
    if args.fourth:
        candidates=[('audit_tolA125',[('reltol=1e-6 abstol=1e-9 vntol=1e-08','reltol=1e-4 abstol=1e-6 vntol=1e-5')]),
                    ('audit_tolB125',[('reltol=1e-6 abstol=1e-9 vntol=1e-08','reltol=1e-5 abstol=1e-6 vntol=1e-6')])]
    if args.shunt:
        candidates=[('audit_shunt1M125',[('rshunt=1e8','rshunt=1e6')]),
                    ('audit_shunt10M125',[('rshunt=1e8','rshunt=1e7')])]
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
        results=list(pool.map(lambda x:candidate(*x),candidates))
    output='data/spice_full_audit_shunt.json' if args.shunt else 'data/spice_full_audit_fourth.json' if args.fourth else 'data/spice_full_audit_third.json' if args.third else 'data/spice_full_audit_second.json' if args.second else 'data/spice_full_audit.json'
    (ROOT/output).write_text(json.dumps(results,indent=2))
