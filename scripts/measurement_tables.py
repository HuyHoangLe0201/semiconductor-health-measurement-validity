"""Export editorial summary tables from existing, unchanged scientific results."""
from pathlib import Path
import hashlib
import json
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
source = ROOT/'data/sil_results.json'
s = json.loads(source.read_text(encoding='utf-8'))
selected = ['clean75', 'adc16', 'adc12', 'skew', 'inputdelay', 'actdelay', 'actdelay_recorded', 'offset2K', 'warming_noT', 'warming_oracle']
labels = {
    'clean75': r'Ideal sensing, 75 $^\circ$C case',
    'adc16': r'16-bit ADC + 2 mA noise',
    'adc12': r'12-bit ADC + 2 mA noise',
    'skew': r'16-bit/noise + 0/1/2 $\mu$s skew',
    'inputdelay': r'16-bit/noise + one-period input delay',
    'actdelay': r'Extra 2 $\mu$s gate delay, undeclared',
    'actdelay_recorded': r'Extra 2 $\mu$s gate delay, recorded',
    'offset2K': r'Injected-record junction offset +2 K',
    'warming_noT': r'Healthy warming, no thermal correction',
    'warming_oracle': r'Healthy warming, oracle correction'
}
conditions = {c['name']: c for c in s['conditions']}
rows = []
audit_rows = []
for name in selected:
    c = conditions[name]
    hv, hr = (1000*x for x in c['contrast'])
    truth = c.get('truth', [0.05, 0.0002])
    ev, er = hv-1000*truth[0], hr-1000*truth[1]
    rows.append(f"{labels[name]}&{hv:.3f}&{ev:+.3f}&{hr:.3f}&{er:+.3f}" + r'\\')
    audit_rows.append({'name': name, 'estimate_SI': c['contrast'], 'truth_SI': truth, 'signed_error_SI': [x-y for x,y in zip(c['contrast'],truth)]})
tex = r"""% Generated from unchanged data/sil_results.json by measurement_tables.py.
\begin{table}[t]\centering\small
\caption{Selected original compiled closed-loop SIL interventions and signed errors. Injected targets are 50 mV and 0.2 m$\Omega$, except for healthy warming (both zero). References are same-condition records except for the explicit warming. Noisy healthy/injected pairs use common random numbers. Oracle compensation uses simulated junction truth.}\label{tab:measurement_sil}
\begin{tabularx}{\textwidth}{@{}Yrrrr@{}}\toprule
Condition&$\widehat h_v$&Error&$\widehat h_r$&Error\\
&\multicolumn{2}{c}{mV}&\multicolumn{2}{c}{m$\Omega$}\\\midrule
""" + '\n'.join(rows) + '\n' + r'\bottomrule\end{tabularx}\end{table}' + '\n'
(ROOT/'measurement_sil_table.tex').write_text(tex, encoding='utf-8')
zsum = norm.ppf(.975) + norm.ppf(.8)
audit = {
    'existing_results_SHA256': hashlib.sha256(source.read_bytes()).hexdigest(),
    'rows': audit_rows,
    'approximate_80_percent_power_standard_error_requirements': {
        'voltage_mV': 50/zsum,
        'resistance_mOhm': .2/zsum,
        'two_sided_alpha': .05,
        'power': .8,
        'scope': 'Per-window Gaussian reference, positive change; approximate two-sided power inversion'
    },
    'thermal_relative_independent_offset_allocations_K': {'voltage': 10/4, 'resistance': .05/.16},
    'SIL_static_anchor_threshold_gain_mV_per_K': 2.5,
    'SIL_static_anchor_10mV_relative_offset_allocation_K': 10/2.5,
    'SIL_resistance_temperature_offset_gain_scope': 'Zero only because the assumed current slopes are temperature independent',
    'ADC_LSB_A': {'16bit': 64/2**16, '12bit': 64/2**12},
    'scope': 'Editorial summaries and conditional allocations only; no new experimental evidence'
}
(ROOT/'validation/measurement_summary_checks.json').write_text(json.dumps(audit,indent=2)+'\n',encoding='utf-8')
print(json.dumps({k:v for k,v in audit.items() if k!='rows'},indent=2))
