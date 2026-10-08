"""Generate manuscript summaries directly from the independent-record audit."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
r=json.loads((ROOT/'data/independent_sil_results.json').read_text())
assert r['calibration_pairs']==200 and r['evaluation_triplets']==400 and r['periods']==4096
scenarios=r['scenarios']
def pc(x): return 100*x['rate']
def interval(x): return f"{pc(x):.2f} [{100*x['exact95'][0]:.2f}, {100*x['exact95'][1]:.2f}]"
v16,r16=scenarios['16']['targets']
macros={'MCVoltageCoverageSixteen':pc(v16['nominal_coverage']),
        'MCResistanceCoverageSixteen':pc(r16['nominal_coverage']),
        'MCVoltageFalseAlarmSixteen':pc(v16['calibrated_false_alarm']),
        'MCResistanceFalseAlarmSixteen':pc(r16['calibrated_false_alarm']),
        'MCResistanceDetectionSixteen':pc(r16['calibrated_detection'])}
(ROOT/'measurement_mc_values.tex').write_text('% Generated from the independent-record audit.\n'+
    '\n'.join('\\newcommand{\\'+key+'}{'+f'{value:.2f}'+'}' for key,value in macros.items())+'\n')
rows=[]
for bits in (16,12):
    a=scenarios[str(bits)]
    assert a['accepted']==400 and a['rejected']==0
    for j,t in enumerate(a['targets']):
        rows.append((bits,j,a,t))
table=[r'\begin{table}[t]\centering\small',
       r'\caption{Held-out independent-record SIL assessment. Each ADC setting uses 200 disjoint null calibration pairs and 400 evaluation triplets. Coverage and decision rates are percentages; nominal coverage is for the injected contrast. Null-calibrated false-alarm and detection use separate targets and a fixed rank threshold. The model uses oracle junction compensation.}\label{tab:independentmc}',
       r'\begin{tabularx}{\textwidth}{@{}Yrrrrr@{}}\toprule',
       r'ADC / target&RMSE&Nominal&Nominal&Calibrated&Calibrated\\',
       r'& &coverage&false alarm&false alarm&detection\\\midrule']
for bits,j,a,t in rows:
    target='voltage' if j==0 else 'resistance'
    table.append(f"{bits}-bit / {target}&{1000*a['RMSE'][j]:.3f}&{pc(t['nominal_coverage']):.2f}&{pc(t['nominal_false_alarm']):.2f}&{pc(t['calibrated_false_alarm']):.2f}&{pc(t['calibrated_detection']):.2f}"+r'\\')
table.extend([r'\bottomrule\end{tabularx}',r'\par\smallskip\footnotesize RMSE units: mV for voltage and m$\Omega$ for resistance. Exact binomial intervals and alternative-coverage diagnostics are in Supplementary Section S5.2.',r'\end{table}'])
(ROOT/'measurement_independent_table.tex').write_text('\n'.join(table)+'\n')

fig,axes=plt.subplots(1,2,figsize=(7.2,3.1),layout='constrained')
for j,ax in enumerate(axes):
    x=np.arange(2)
    for shift,key,label,mark in [(-.14,'nominal_false_alarm','Nominal false alarm','o'),
                                (0,'calibrated_false_alarm','Null-calibrated false alarm','s'),
                                (.14,'calibrated_detection','Null-calibrated detection','^')]:
        v=[scenarios[str(bits)]['targets'][j][key] for bits in (16,12)]
        y=np.array([pc(z) for z in v]);lo=np.array([100*z['exact95'][0] for z in v]);hi=np.array([100*z['exact95'][1] for z in v])
        ax.errorbar(x+shift,y,yerr=np.array([y-lo,hi-y]),fmt=mark,markersize=4,capsize=3,label=label)
    ax.axhline(5,color='gray',ls='--',lw=.8);ax.set_xticks(x,['16 bit','12 bit'])
    ax.set(title='50 mV voltage change' if j==0 else r'0.2 m$\Omega$ resistance change',ylabel='Held-out rate (%)',ylim=(-3,103))
    ax.grid(alpha=.2)
axes[0].legend(fontsize=7,loc='center right')
fig.savefig(ROOT/'figures/independent_sil_decision.pdf',bbox_inches='tight')
fig.savefig(ROOT/'figures/independent_sil_decision.png',dpi=180,bbox_inches='tight');plt.close(fig)

methods=r'''\subsection{Independent-record acquisition and held-out decisions}\label{sec:independentmc}
The additional experiment fixes one acquisition context throughout: a T-type drive with 300 V DC link, 2 mH inductance, 0.25 $\Omega$ resistance, 50 $\mu$s period, 1 $\mu$s blanking and a 75 $^\circ$C case boundary. Records span 4096 periods (204.8 ms). The modeled health coefficients refer to $T_0=25$ $^\circ$C; their injected increments are assumed temperature independent. The retained measured-current endpoints exceed 0.5 A in magnitude and the plant is constrained below 31.99 A, within the 64 A ADC span. This is the declared local surrogate range, not an empirically calibrated current--voltage domain. Current noise is 2 mA per channel; ADC resolution is either 16 or 12 bits. Each record regenerates its feedback, gates and electrothermal plant. Machine scale is supplied exactly and junction compensation remains an oracle. The intended changes are 50 mV and 0.2 m$\Omega$; a second healthy record defines a genuine zero-change comparison.

For each ADC setting, 200 healthy/healthy pairs form a calibration partition, followed by 400 evaluation triplets consisting of a healthy reference, a separate healthy test and an injected test. All three current-noise streams are independent. The reference is shared only between the two comparisons within an evaluation triplet; trials and calibration/evaluation partitions are disjoint. The nominal target covariance includes the two record contributions in~\eqref{eq:baselinecov}. Its midpoint MA(1) weight uses $\sigma_{\mathrm{wt}}^2=(2\ \mathrm{mA})^2+\mathrm{LSB}^2/12$. The additive uniform-quantization variance is a reference approximation, not an exact likelihood for quantized feedback. Neither label nor validation error enters the record fit.

For each target separately, calibration scores are $S_i=|\widehat{\Delta h}_i|/s_{\Delta,i}$ under the null. The decision rejects zero change when $|\widehat{\Delta h}|>q s_\Delta$, where $q$ is the 191st ordered score among the 200 calibration pairs. This is the finite-rank calibration rule $\lceil(n_{\mathrm{cal}}+1)(1-\alpha)\rceil$ for $\alpha=0.05$ \cite{angelopoulos2022}. Under exchangeable null pairs it bounds the marginal false-alarm probability over calibration and a fresh null pair; it does not guarantee 5\% conditional error for every realized calibration partition. It supplies no alternative-coverage guarantee, simultaneous-target guarantee or protection against operating-condition shifts. Injected-record coverage and detection are therefore measured on the held-out partition.

\input{measurement_independent_table}
\begin{figure}[t]\centering
\includegraphics[width=\textwidth]{independent_sil_decision.pdf}
\caption{Independent-record closed-loop decisions. Points and exact 95\% binomial intervals use 400 held-out null or injected comparisons per ADC setting; the dashed line marks 5\%. Calibration uses 200 different null pairs. The circuit, thermal compensation and sensing remain software models.\AIFigureNote}\label{fig:independentdecision}
\end{figure}
'''
results=[]
for bits in (16,12):
    a=scenarios[str(bits)];v,t=a['targets']
    results.append(f"For {bits}-bit acquisition, nominal injected-target coverage is {pc(v['nominal_coverage']):.2f}\\% for voltage and {pc(t['nominal_coverage']):.2f}\\% for resistance. The held-out null-calibrated false-alarm rates are {pc(v['calibrated_false_alarm']):.2f}\\% and {pc(t['calibrated_false_alarm']):.2f}\\%, and detection rates are {pc(v['calibrated_detection']):.2f}\\% and {pc(t['calibrated_detection']):.2f}\\%, respectively.")
comparison=[]
for bits in (16,12):
    a=scenarios[str(bits)]['common_noise_comparator']
    comparison.append(f"{bits} bits: {a['SD_ratio'][0]:.3f} for voltage and {a['SD_ratio'][1]:.3f} for resistance")
requirement=f"The 16-bit median independent-difference standard errors are {1000*scenarios['16']['median_nominal_SE'][0]:.3f} mV and {1000*scenarios['16']['median_nominal_SE'][1]:.3f} m$\\Omega$. The voltage precision is below the illustrative 17.85 mV screening allocation, whereas the resistance precision exceeds 0.0714 m$\\Omega$. The held-out detection rates assess the actual feedback decision; the Gaussian screening values alone are not substituted for that assessment."
end=' '.join(results)+r''' All 400 evaluation triplets per setting are accepted. These rates assess this fixed software acquisition context; they are not achieved hardware specifications. Figure~\ref{fig:independentdecision} connects a defined record difference, uncertainty weight, calibration partition and held-out decision in one model, rather than transferring the abstract-design Gaussian power to feedback.
'''+ '\n\n'+requirement+r'''

The first 100 evaluation references also have an injected counterpart driven by their same current-noise stream, making a direct common-versus-independent pairing comparison. The common/independent standard-deviation ratios are '''+'; '.join(comparison)+r'''. Pairing can change the dispersion and must be declared. Supplementary Section S5.2 reports bias, standard errors, alternative coverage, exact rate intervals, parity checks, seeds and source hashes. No null calibration removes an unmodeled machine-scale change, temperature-reference drift or sampling error.
'''
(ROOT/'measurement_independent_sil.tex').write_text(methods+'\n'+end)

supp=r'''\subsection{Independent-record Monte Carlo and calibration partition}\label{sec:independentdetails}
This extension uses the circuit, controller, path operator and 36-dimensional quotient of Section S5.1 without changing the original source units. A separate native batching harness reproduces synchronous acquisition, the four-substep integration grid, NumPy nearest-even ADC rounding and oracle junction compensation. Two full 4096-period parity records, at 16 and 12 bits, are compared with the original Python orchestration. Current, indication and junction differences are below $10^{-6}$ in their respective units and all selected gates coincide. Python SVD independently checks the complete-operator rank, estimates and full $2\times2$ target covariance against the new QR covariance extraction; numerical differences are stored in the audit.

The experiment has 200 null calibration pairs and 400 evaluation triplets per ADC setting. Each trial starts from the same cold-modal state at a fixed 75 $^\circ$C case boundary and regenerates feedback for 204.8 ms. Reference, null-test and injected-test records use independent Gaussian current-noise arrays, each with 2 mA per-channel standard deviation. Noise, gate and current-array hashes are recorded. PCG64 seeds are $202610080000+10^6b+10i+j$, where $b\in\{16,12\}$ is ADC resolution, $i=0,\ldots,599$ is trial index, and $j=0,1,2$ denotes reference, null and injected records. Indices 0--199 are calibration; 200--599 are evaluation. Calibration has no injected record. The first 100 evaluation trials also run an injected record with $j=0$ solely for the explicitly labeled common-random-number comparison. Noise independence applies across records with distinct seeds, not to the two contrasts that share one reference inside a triplet.

The nominal covariance sums the two fitted record covariance matrices, with no random-error cross term because the primary streams are independent. Both estimates and nominal uncertainty use the reconstructed measured regressors and selected-period covariance. The weight has variance $(2\ \mathrm{mA})^2+(64/2^b)^2/12$; quantization independence and the MA(1) law are approximations used for weighting. Feedback/regressor dependence is not silently converted into a valid Gaussian likelihood. Rank-deficient or out-of-range records abort or are counted as rejected; no such record is discarded to improve rates. In the completed evaluation no rejection occurs.

The 191st absolute standardized null score defines a separate decision threshold for each ADC/target combination. The partition and rank rule are specified in the run script before acquisition. The thresholds are computed using only the designated calibration records and applied without retuning to the held-out records; numerical thresholds are calculated after acquisition, not used to choose the generated trajectories. This is an established finite-rank exchangeability construction \cite{angelopoulos2022}; it is used here as a transparent decision assessment, not presented as a new method. The guarantee is marginal over simulated calibration and null pairs under the same acquisition distribution. It does not extend to injected-target coverage or repeated overlapping online windows. The nominal 1.96-standard-error interval and the null-calibrated interval are both evaluated empirically around the injected truth. Exact Clopper--Pearson intervals accompany every reported rate, conditional on the realized calibration partition; they do not include the variability of refitting the calibration threshold. All 400 null comparisons are independent across trials, as are all 400 injected comparisons; those two collections share references pairwise and are not combined as 800 independent observations.
'''
detail=[r'\begin{table}[t]\centering\small',r'\caption{Independent-record SIL error and uncertainty diagnostics. Units are mV for voltage and m$\Omega$ for resistance; half-width uses the frozen null-calibrated multiplier.}\label{tab:independenterrors}',r'\begin{tabular}{@{}llrrrrr@{}}\toprule ADC&Target&Bias&SD&RMSE&Median SE&Half-width\\\midrule']
for bits,j,a,t in rows:
    detail.append(f"{bits} bit&{'Voltage' if j==0 else 'Resistance'}&"+'&'.join(f'{1000*a[key][j]:.3f}' for key in ['bias','SD','RMSE','median_nominal_SE','median_calibrated_halfwidth'])+r'\\')
detail += [r'\bottomrule\end{tabular}\end{table}']
for key,label in [('nominal_coverage','Nominal injected coverage'),('calibrated_alternative_coverage','Null-calibrated injected coverage'),('nominal_false_alarm','Nominal false alarm'),('calibrated_false_alarm','Null-calibrated false alarm'),('nominal_detection','Nominal detection'),('calibrated_detection','Null-calibrated detection')]:
    vals=[f"{bits}-bit {'voltage' if j==0 else 'resistance'}: {interval(t[key])}\\%" for bits,j,a,t in rows]
    detail.append(label+' (rate [exact 95\\% interval]): '+ '; '.join(vals)+'.\n')
detail.append('Frozen calibration multipliers (voltage, resistance): '+ '; '.join(f"{bits} bits ({scenarios[str(bits)]['calibration_quantiles'][0]:.4f}, {scenarios[str(bits)]['calibration_quantiles'][1]:.4f})" for bits in (16,12))+'.\n')
detail.append(r'''\begin{figure}[t]\centering\includegraphics[width=\textwidth]{independent_sil_decision.pdf}
\caption{Held-out decision rates in the independent-record SIL study. Each rate has 400 evaluation comparisons and an exact binomial interval. The thresholds are calibrated using 200 disjoint null pairs per ADC setting.\AIFigureNote}\label{fig:independentdecision}\end{figure}
The reproducibility package supplies every trial's seed, record estimates, full target covariance, status and noise/current/gate hashes in \texttt{data/independent\_sil/trials.jsonl}. Complete arrays are retained for parity records and representative calibration/evaluation trials. Remaining traces are reproducible from the seeds and native sources; they are not claimed to be all stored. The original common-noise evidence is preserved separately. All execution remains offline on a PC, with ideal machine references and oracle junction compensation. No physical calibration, HIL or converter-bench data are introduced.
''')
(ROOT/'measurement_independent_supplement.tex').write_text(supp+'\n'+'\n'.join(detail))
print('Independent-record tables, figures and manuscript sections generated from completed audit.')
