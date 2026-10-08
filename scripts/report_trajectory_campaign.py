from pathlib import Path
import json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'data/trajectory_campaign/results.json').read_text());recs=d['records'];windows=d['contract']['windows']
plt.rcParams.update({'font.family':'DejaVu Serif','font.size':9,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
fig,axs=plt.subplots(2,2,figsize=(10,7.0),layout='constrained')
ax=axs[0,0]
for j,marker,label in [(0,'o',r'$h_v$'),(1,'s',r'$h_r$')]:
 vals=[100*d['summary'][str(n)]['support_counts'][j]/144 for n in windows]
 ax.plot(windows,vals,color='black',ls=('-','--')[j],marker=marker,label=label)
 for n,v in zip(windows,vals):
  offset=(8 if j==0 else -14) if n==128 else (-14,7)[j]
  ax.annotate(f'{v:.1f}',(n,v),xytext=(0,offset),textcoords='offset points',ha='center',fontsize=8)
ax.set(xscale='log',xticks=windows,xticklabels=windows,ylim=(0,110),xlabel='Window (control periods)',ylabel='Supported cases (%)',title='(a) Test each target before assigning precision')
ax.legend(loc='lower right');ax.grid(axis='y',color='.85')
ax=axs[0,1];ranks=np.array([[w['rank'] for w in r['windows']] for r in recs])
im=ax.imshow(ranks,aspect='auto',cmap='Greys',vmin=0,vmax=36,interpolation='nearest')
ax.set(xticks=range(4),xticklabels=windows,yticks=[0,23,47,71,95,119,143],yticklabels=[1,24,48,72,96,120,144],xlabel='Window (control periods)',ylabel='Predeclared case index',title='(b) Full 50-column operator; maximum rank 36')
fig.colorbar(im,ax=ax,label='Retained numerical rank',ticks=[0,12,24,36],shrink=.85)
ax=axs[1,0]
for j,marker,label in [(0,'o',r'$h_v$/50 mV'),(1,'s',r'$h_r$/0.2 m$\Omega$')]:
 for a,A in enumerate((2,8,20)):
  vals=[r['windows'][-1]['targets'][j]['se']/(.05,.0002)[j] for r in recs if r['config']['A']==A and r['windows'][-1]['targets'][j]['supported']]
  xpos=a+(-.14,.14)[j];ax.scatter(xpos+np.linspace(-.075,.075,len(vals)),vals,s=12,marker=marker,facecolors=('none','black')[j],edgecolors='black',linewidths=.6,label=label if a==0 else None)
  ax.plot([xpos-.09,xpos+.09],[np.median(vals)]*2,color='black',lw=2)
ax.axhline(1,color='.4',ls=':',lw=1)
ax.set(yscale='log',xticks=range(3),xticklabels=[2,8,20],xlabel='Nominal reference amplitude (A)',ylabel='Reference standard error / target change',title='(c) Supported does not imply useful precision')
ax.legend(fontsize=8,loc='upper right');ax.grid(axis='y',color='.9')
ax=axs[1,1]
for j,marker,label in [(0,'o',r'$h_v$'),(1,'+',r'$h_r$')]:
 ws=[r['windows'][-1] for r in recs if r['windows'][-1]['targets'][j]['supported']]
 ax.scatter([w['smin'] for w in ws],[w['targets'][j]['kappa'] for w in ws],s=17,marker=marker,facecolors='none' if j==0 else 'black',edgecolors='black',linewidths=.7,label=label)
ax.set(xscale='log',yscale='log',xlabel=r'Smallest retained $\sigma_+(XS/\sqrt{m})$',ylabel=r'Target condition $\kappa_l$',title='(d) Conditioning belongs to the target and metric')
ax.legend(loc='upper right');ax.grid(color='.9')
fig.savefig(ROOT/'figures/trajectory_support.pdf',bbox_inches='tight');fig.savefig(ROOT/'figures/trajectory_support.png',dpi=180,bbox_inches='tight')
print('Figure saved; unsupported cases omitted from precision panels only.')
