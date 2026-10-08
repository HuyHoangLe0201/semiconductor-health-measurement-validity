"""Reanalysis of existing NASA PCoE device files; no converter hardware claim.

Original files are kept intact. Temperatures are package channels, not junction
truth. All processing rules are fixed without fitting to a held-out waveform.
"""
from pathlib import Path
import argparse,hashlib,json,shutil,gc
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'external'/'nasa';DATA=ROOT/'data';FIG=ROOT/'figures'
SOURCE_URL='https://data.nasa.gov/dataset/insulated-gate-bipolar-transistor-igbt-accelerated-aging'

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def extract(source):
    RAW.mkdir(parents=True,exist_ok=True)
    manifest=[]; rows=[]; bins=[];exclusions=[]
    for device in range(2,6):
        name=f'Device{device}__1.mat';src=source/name;target=RAW/name
        if not target.exists():shutil.copy2(src,target)
        assert digest(src)==digest(target)
        m=loadmat(target,simplify_cells=True)['measurement']
        manifest.append(dict(device=device,file=name,bytes=target.stat().st_size,
                             sha256=digest(target),local_source=str(src),
                             original_remote_archive_verified=False,total_transient_captures=len(m['transient'])))
        steady=m['steadyState'];st=np.array([x['timeEpoch'] for x in steady])
        package=np.array([x['timeDomain']['packageTemperature'] for x in steady],float)
        good=(st>0)&np.isfinite(package)
        st=st[good];package=package[good];order=np.argsort(st);st=st[order];package=package[order]
        for cycle,x in enumerate(m['transient']):
            td=x['timeDomain'];dt=float(td['dt']);t=float(x['timeEpoch'])
            i=np.asarray(td['collectorEmitterCurrentSignal'],float)
            v=np.asarray(td['collectorEmitterVoltage'],float)
            gate=np.asarray(td['gateEmitterVoltage'],float)
            n=len(i);assert len(v)==n==len(gate)
            on=gate>10.;edges=np.flatnonzero(np.diff(on.astype(int)))+1
            # Exclude 10 us around every gate boundary and capture endpoints.
            guard=int(np.ceil(10e-6/dt));valid=on&(i>1)&(v>0)&np.isfinite(i+v)
            valid[:guard]=False;valid[-guard:]=False
            for edge in edges:valid[max(0,edge-guard):min(n,edge+guard)]=False
            # Use nonoverlapping contiguous 2-us medians, never individual samples
            # as independent replicates. Preserve waveform and block identifiers.
            width=int(round(2e-6/dt));bi=[];bv=[]
            for start in range(0,n-width+1,width):
                if np.all(valid[start:start+width]):
                    im=float(np.median(i[start:start+width]));vm=float(np.median(v[start:start+width]))
                    bi.append(im);bv.append(vm)
                    bins.append((device,cycle+1,start*dt,im,vm))
            bi=np.array(bi);bv=np.array(bv)
            idx=int(np.argmin(abs(st-t)));distance=float(abs(st[idx]-t))
            temp=float(package[idx]) if distance<=2. else np.nan
            if len(bi)<10:
                exclusions.append(dict(device=device,capture=cycle+1,reason='fewer than 10 eligible 2-us blocks',eligible_blocks=len(bi)))
                continue
            X=np.c_[np.ones(len(bi)),bi];b=np.linalg.lstsq(X,bv,rcond=None)[0]
            rmse=float(np.sqrt(np.mean((bv-X@b)**2)))
            rows.append((device,cycle+1,t,len(bi),float(np.median(bi)),float(np.median(bv)),
                         float(np.quantile(bi,.05)),float(np.quantile(bi,.95)),float(b[0]),float(b[1]),rmse,temp,distance))
        print('NASA_EXTRACTED',device,len(m['transient']),flush=True)
        del m;gc.collect()
    rows=np.array(rows);bins=np.array(bins)
    np.savetxt(DATA/'nasa_waveform_summary.csv',rows,delimiter=',',header='device,capture,source_epoch,blocks,median_current_channel,median_Vce_V,I05,I95,affine_intercept_V,affine_slope_per_current_unit,within_waveform_rmse_V,nearest_package_temperature_C,temperature_time_distance_s',comments='')
    np.savez_compressed(DATA/'nasa_onstate_blocks.npz',blocks=bins,summary=rows)
    (RAW/'provenance.json').write_text(json.dumps(dict(source_url=SOURCE_URL,evidence_kind='third_party_measured_device_data',files=manifest,
       current_channel='collectorEmitterCurrentSignal, preserved as stored; SI calibration not independently established',temperature_channel='packageTemperature; nearest record <=2 s; not measured junction',
       excluded=exclusions,procedure='gateEmitterVoltage >10 V; positive current channel >1; Vce>0; exclude 10 us at edges and endpoints; nonoverlapping 2 us medians; at least 10 blocks per waveform'),indent=2),encoding='utf-8')
    return rows,bins

def report(rows,bins):
    summary=[];diagnostics=[]
    for a in rows:
        b=bins[(bins[:,0]==a[0])&(bins[:,1]==a[1])]
        current,voltage=b[:,3],b[:,4]
        cut=int(.6*len(b));train=np.arange(cut);test=np.arange(cut+2,len(b))
        if len(test)<4:continue
        anchor=current[train].mean()
        X=np.c_[np.ones(len(b)),current-anchor]
        coef=np.linalg.lstsq(X[train],voltage[train],rcond=None)[0]
        fitall=np.linalg.lstsq(X,voltage,rcond=None)[0]
        s=np.linalg.svd(np.c_[np.ones(len(b)),current]/np.linalg.norm(np.c_[np.ones(len(b)),current],axis=0),compute_uv=False)
        condition=float(s[0]/s[-1])
        affine=float(np.sqrt(np.mean((voltage[test]-X[test]@coef)**2)))
        constant=float(np.sqrt(np.mean((voltage[test]-voltage[train].mean())**2)))
        diagnostics.append(dict(device=int(a[0]),capture=int(a[1]),blocks=len(b),
          condition_column_normalized=condition,current_anchor=anchor,heldout_affine_rmse_V=affine,
          heldout_constant_rmse_V=constant,local_drop_fit_change_V=float(abs(coef[0]-fitall[0])),
          zero_current_intercept_fit_change_V=float(abs((coef[0]-anchor*coef[1])-(fitall[0]-anchor*fitall[1])))))
    for device in range(2,6):
        a=rows[rows[:,0]==device];b=bins[bins[:,0]==device]
        ds=[x for x in diagnostics if x['device']==device]
        # This is an out-of-time discrepancy check, not a causal aging model.
        split=int(np.ceil(.5*len(a)));train=a[:split];test=a[split:]
        temporal={}
        for label,cols in [('current_only',[4]),('current_and_package_temperature',[4,11])]:
            center=train[:,cols].mean(axis=0);scale=train[:,cols].std(axis=0)
            Xtr=np.c_[np.ones(len(train)),(train[:,cols]-center)/scale]
            Xte=np.c_[np.ones(len(test)),(test[:,cols]-center)/scale]
            coef=np.linalg.lstsq(Xtr,train[:,5],rcond=None)[0]
            temporal[label]=dict(train_captures=len(train),test_captures=len(test),
                heldout_rmse_V=float(np.sqrt(np.mean((Xte@coef-test[:,5])**2))),
                train_design_condition=float(np.linalg.cond(Xtr)))
        summary.append(dict(device=device,captures=len(a),blocks=len(b),
           package_temp_range_C=[float(np.nanmin(a[:,11])),float(np.nanmax(a[:,11]))],
           matched_temp_captures=int(np.isfinite(a[:,11]).sum()),
           current_channel_range=[float(np.min(b[:,3])),float(np.max(b[:,3]))],
           median_current_span_90=float(np.median(a[:,7]-a[:,6])),
           median_affine_rmse_mV=float(np.median(a[:,10])*1000),
           intercept_range_V=[float(a[:,8].min()),float(a[:,8].max())],
           slope_range=[float(a[:,9].min()),float(a[:,9].max())],
           split_diagnostic_captures=len(ds),median_condition_column_normalized=float(np.median([x['condition_column_normalized'] for x in ds])),
           median_heldout_affine_rmse_mV=float(np.median([x['heldout_affine_rmse_V'] for x in ds])*1000),
           median_heldout_constant_rmse_mV=float(np.median([x['heldout_constant_rmse_V'] for x in ds])*1000),
           median_intercept_split_sensitivity_mV=float(np.median([x['zero_current_intercept_fit_change_V'] for x in ds])*1000),
           median_local_drop_split_sensitivity_mV=float(np.median([x['local_drop_fit_change_V'] for x in ds])*1000),
           capture_time_temperature_distance_max_s=float(a[:,12].max()),
           voltage_range_V=[float(b[:,4].min()),float(b[:,4].max())],
           chronological_capture_prediction=temporal))
    print(json.dumps(summary,indent=2),flush=True)
    result=dict(evidence_kind='third_party_measured_device_data',hardware_run_by_authors=False,
       upstream=SOURCE_URL,devices=summary,diagnostics=diagnostics,
       total_eligible_captures=len(rows),total_blocks=len(bins),total_split_diagnostic_captures=len(diagnostics),
       notes='No healthy labels inferred from chronology; no causal temperature coefficient or converter accuracy claimed. Current channel calibration not independently verified. Package temperature is a nearest-time contextual channel, not junction truth. Held-out blocks are chronological and separated by two blocks; block dependence remains, so no iid confidence intervals are used.')
    (DATA/'external_device_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(2,2,figsize=(8.6,6.1),constrained_layout=True)
    for ax,d in zip(axes.flat,range(2,6)):
        a=rows[rows[:,0]==d]
        im=ax.scatter(a[:,1],a[:,5],c=a[:,11],cmap='plasma',vmin=25,vmax=285,s=15)
        ax.set(xlabel='Source capture index',ylabel='Median gated Vce (V)',title=f'NASA device {d}')
        ax.grid(alpha=.2)
    fig.colorbar(im,ax=axes.ravel().tolist(),label='Nearest package temperature (deg C)',shrink=.83)
    fig.savefig(FIG/'external_device_trace.pdf');fig.savefig(FIG/'external_device_trace.png',dpi=200);plt.close(fig)
    table=['\\begin{table}[t]\\centering\\small','\\caption{Reanalysis of third-party NASA device measurements. Current span uses the stored channel units; voltage errors are mV. Condition numbers use column-normalized $[1,I]$. Held-out errors use the last 40\\% of each gated waveform. No row establishes converter-level health accuracy.}\\label{tab:nasa}',
      '\\begin{tabular}{@{}lrrrrr@{}}\\toprule','Device&Captures&Median span&Condition&Affine RMSE&Constant RMSE\\\\\\midrule']
    for s in summary:
        table.append(f"{s['device']}&{s['captures']}&{s['median_current_span_90']:.3f}&{s['median_condition_column_normalized']:.0f}&{s['median_heldout_affine_rmse_mV']:.2f}&{s['median_heldout_constant_rmse_mV']:.2f}\\\\")
    table += ['\\bottomrule\\end{tabular}','\\end{table}']
    (FIG/'external_device_summary.tex').write_text('\n'.join(table)+'\n',encoding='ascii')
    stats='The selection retains %d of 423 source captures and %d median blocks. Median 5th--95th percentile current spans are %.3f--%.3f stored units, with column-normalized design condition numbers %.0f--%.0f. Across devices, median held-out affine errors are %.2f--%.2f mV, whereas the constant-voltage reference gives %.2f--%.2f mV. Thus a small local prediction error does not establish a benefit from fitting both coefficients. Changing from the leading-60\\%% fit to the complete waveform changes the zero-current intercept by a device median of %.2f--%.2f mV, compared with %.2f--%.2f mV for voltage evaluated at the training current mean. These are split sensitivities, not confidence intervals or aging magnitudes.\n' % (
        len(rows),len(bins),min(s['median_current_span_90'] for s in summary),max(s['median_current_span_90'] for s in summary),
        min(s['median_condition_column_normalized'] for s in summary),max(s['median_condition_column_normalized'] for s in summary),
        min(s['median_heldout_affine_rmse_mV'] for s in summary),max(s['median_heldout_affine_rmse_mV'] for s in summary),
        min(s['median_heldout_constant_rmse_mV'] for s in summary),max(s['median_heldout_constant_rmse_mV'] for s in summary),
        min(s['median_intercept_split_sensitivity_mV'] for s in summary),max(s['median_intercept_split_sensitivity_mV'] for s in summary),
        min(s['median_local_drop_split_sensitivity_mV'] for s in summary),max(s['median_local_drop_split_sensitivity_mV'] for s in summary))
    stats += 'One ten-block capture has fewer than four withheld blocks after the split guard and is excluded from the split diagnostics; the other 309 captures contribute to those summaries.\n'
    (ROOT/'nasa_computed_summary.tex').write_text(stats,encoding='ascii')

if __name__=='__main__':
    default_source=RAW
    parser=argparse.ArgumentParser();parser.add_argument('--source',type=Path,default=default_source);parser.add_argument('--reuse',action='store_true');args=parser.parse_args()
    if args.reuse:
        z=np.load(DATA/'nasa_onstate_blocks.npz');rows,bins=z['summary'],z['blocks']
    else:rows,bins=extract(args.source)
    report(rows,bins)
