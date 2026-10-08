"""Optional graphical abstract: a vector schematic, not new experimental evidence."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'data/independent_sil_results.json').read_text())
decision=data['scenarios']['16']['targets']
plt.rcParams.update({'font.family':'DejaVu Sans','svg.fonttype':'none'})
plt.rcParams['pdf.fonttype']=42
fig,ax=plt.subplots(figsize=(15,6),dpi=300)
fig.patch.set_facecolor('white')
ax.set(xlim=(0,15),ylim=(0,5.4));ax.axis('off')
colors=['#eef5fb','#eef5fb','#eaf5f1','#eaf5f1']
items=[
    ('Acquire', 'Currents and gates\nTiming / reference'),
    ('Model', 'Conduction paths\nShared influences'),
    ('Select', 'Health contrasts\nThermal baseline'),
    ('Report', 'Uncertainty\nDetection power')
]
xs=[.15,3.92,7.69,11.46]
for x,color,(title,body) in zip(xs,colors,items):
    box=FancyBboxPatch((x,2.8),3.38,2.02,boxstyle='round,pad=0.05,rounding_size=0.12',facecolor=color,edgecolor='#385268',linewidth=1.4)
    ax.add_patch(box)
    ax.text(x+1.69,4.38,title,ha='center',va='center',fontsize=23,fontweight='bold',color='#243d50')
    ax.text(x+1.69,3.57,body,ha='center',va='center',fontsize=20,linespacing=1.65,color='#243d50')
for a,b in zip(xs[:-1],xs[1:]):
    ax.add_patch(FancyArrowPatch((a+3.45,3.85),(b-.10,3.85),arrowstyle='-|>',mutation_scale=15,linewidth=1.4,color='#385268'))
ax.text(7.5,2.26,'Structural support does not establish measurement accuracy',ha='center',fontsize=21,fontweight='bold',color='#243d50')
voltage=100*decision[0]['calibrated_detection']['rate']
resistance=100*decision[1]['calibrated_detection']['rate']
ax.text(3.7,1.52,'Independent-record software study\n200 null calibration pairs\n400 held-out triplets per ADC',ha='center',va='center',fontsize=17,linespacing=1.65,color='#385268')
ax.text(11.3,1.52,f'16-bit observed detection\n50 mV: {voltage:.2f}%\n0.2 mOhm: {resistance:.2f}%',ha='center',va='center',fontsize=17,linespacing=1.65,color='#385268')
ax.text(7.5,.33,'Simulated acquisition with oracle thermal compensation; physical accuracy unverified.',ha='center',fontsize=14,color='#52616b')
ax.text(7.5,.04,'Visualization code prepared with assistance from OpenAI Codex.',ha='center',fontsize=11,color='#52616b')
fig.subplots_adjust(left=.01,right=.99,top=.99,bottom=.01)
fig.savefig(ROOT/'figures/measurement_graphical_abstract.png',dpi=300,facecolor='white')
fig.savefig(ROOT/'figures/measurement_graphical_abstract.svg',facecolor='white')
fig.savefig(ROOT/'figures/measurement_graphical_abstract.pdf',facecolor='white',metadata={'Author':'','Title':'Graphical abstract: identifiable health contrasts and measurement uncertainty'})
out=ROOT.parent/'output/pdf';out.mkdir(parents=True,exist_ok=True)
(out/'P1_Measurement_Graphical_Abstract.pdf').write_bytes((ROOT/'figures/measurement_graphical_abstract.pdf').read_bytes())
plt.close(fig)
print('Optional graphical abstract saved as PDF, PNG and SVG.')
