"""Build a configurable derivative, preserving all original SIL sources."""
from pathlib import Path
import hashlib, json, subprocess
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'validation/trajectory_native'
OUT.mkdir(parents=True,exist_ok=True)
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
changes={
 'sil_core.c':[
 ('#include "sil_basis.h"','#include "sil_basis.h"\n#include "tv_config.h"'),
 ('16+5*sin(2*3.141592653589793*7*(t+50e-6))+3*sin(2*3.141592653589793*13*(t+50e-6))','tv_amp*(1+tv_mod*(5*sin(2*3.141592653589793*7*(t+50e-6))+3*sin(2*3.141592653589793*13*(t+50e-6)))/16)'),
 ('20*sin(2*3.141592653589793*25','tv_E*sin(2*3.141592653589793*25'),
 ('.025*(u[j]-mean-.25*sense[j])','(50e-6/tv_L)*(u[j]-mean-tv_R*sense[j])'),
 ('+shift[j]-.45','+shift[j]+tv_phase'),
 ('-.25*mid[p]','-tv_R*mid[p]'),
 ('40*(hi[p]-lo[p])','(tv_L/50e-6)*(hi[p]-lo[p])')],
 'sil_plant.c':[
 ('#include <string.h>','#include <string.h>\n#include "tv_config.h"'),
 ('-20*sin(157.07963267948966','-tv_E*sin(157.07963267948966'),
 ('-.25*current[phase])/.002','-tv_R*current[phase])/tv_L')]
}
provenance=[]
for source,replacements in changes.items():
 p=ROOT/'scripts'/source; s=p.read_text()
 for a,b in replacements:
  assert a in s,(source,a)
  s=s.replace(a,b)
 dest=OUT/source.replace('sil_','tv_'); dest.write_text(s)
 provenance.append(dict(original=str(p.relative_to(ROOT)),original_sha256=digest(p),derived=str(dest.relative_to(ROOT)),derived_sha256=digest(dest),substitutions=replacements))
(OUT/'sil_basis.h').write_bytes((ROOT/'scripts/sil_basis.h').read_bytes())
(OUT/'tv_config.h').write_text('extern double tv_L,tv_R,tv_E,tv_amp,tv_mod,tv_phase;\n')
driver=(ROOT/'scripts/sil_mc_native.c').read_text().split('/* Additional QR fit')[0]
driver=driver.replace('#include "sil_core.c"','#include "tv_core.c"\ndouble tv_L=.002,tv_R=.25,tv_E=20,tv_amp=16,tv_mod=1,tv_phase=-.45;\nAPI void tv_set(double L,double R,double E,double A,double m,double phase){tv_L=L;tv_R=R;tv_E=E;tv_amp=A;tv_mod=m;tv_phase=phase;}')
driver=driver.replace('for(h=0;h<4;h++)sil_advance(k*50e-6+h*.25e-6,state,freeq,.25e-6,Tc,aged,1);','for(h=0;h<tv_refine;h++)sil_advance(k*50e-6+h*1e-6/tv_refine,state,freeq,1e-6/tv_refine,Tc,aged,1);')
driver=driver.replace('for(h=0;h<4;h++)sil_advance(k*50e-6+1e-6+h*12.25e-6,state,q,12.25e-6,Tc,aged,1);','for(h=0;h<tv_refine;h++)sil_advance(k*50e-6+1e-6+h*49e-6/tv_refine,state,q,49e-6/tv_refine,Tc,aged,1);')
driver=driver.replace('API int sil_mc_run','static int tv_refine=4;\nAPI void tv_set_refine(int r){tv_refine=r;}\nAPI int sil_mc_run')
driver=driver.replace('sil_control(k*50e-6,sense,old,q);','if(tv_replay)memcpy(q,states+3*k,3*sizeof(int));else sil_control(k*50e-6,sense,old,q);')
driver=driver.replace('static int tv_refine=4;','static int tv_refine=4,tv_replay=0;\nAPI void tv_set_replay(int r){tv_replay=r;}')
(OUT/'tv_driver.c').write_text(driver)
plant=OUT/'tv_plant.c'; plant.write_text(plant.read_text()+'\nAPI void tv_ode(double t,const double *x,const int *q,double Tc,double *out){ode(t,x,q,Tc,0,1,out);}\n')
for unit in provenance: unit['derived_sha256']=digest(ROOT/unit['derived'])
bat=OUT/'compile.cmd'
bat.write_text('@echo off\ncall "C:\\Program Files (x86)\\Microsoft Visual Studio\\2019\\BuildTools\\VC\\Auxiliary\\Build\\vcvars64.bat" >nul\ncl /nologo /O2 /fp:precise /W4 /LD tv_driver.c tv_plant.c /link /OUT:trajectory.dll\n')
result=subprocess.run(['cmd','/c',str(bat)],cwd=OUT,capture_output=True,text=True)
(OUT/'build.log').write_text(result.stdout+result.stderr)
assert result.returncode==0,result.stdout+result.stderr
(OUT/'source_provenance.json').write_text(json.dumps(dict(original_units_preserved=True,unused_fixed_basis_fit=True,units=provenance,build_returncode=result.returncode),indent=2))
print(OUT/'trajectory.dll')
