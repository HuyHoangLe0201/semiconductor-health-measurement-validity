"""Geometric companion plates. Analytic scenarios, no invented measurements."""
def append_visual_extensions(add, picture, panel):
    b=picture()+panel(.1,15.45,7.6,'(a) Equal sample count and mean current')+panel(8.35,15.45,8.05,'(b) Profile out the slope column')
    b+=r'''\begin{axis}[mathaxis,at={(.8cm,10.5cm)},anchor=south west,width=7.2cm,height=4.5cm,
 xmin=-.15,xmax=2.15,ymin=.1,ymax=1.15,xtick={0,.8,1,1.2,2},ytick={.35,.85},yticklabels={Wide,Narrow},
 xlabel={Normalized current $u$},ylabel={Design},clip=false]
\addplot[only marks,mark=o,mark size=3pt] coordinates {(.8,.85)(1.2,.85)};
\addplot[only marks,mark=*,mark size=2pt] coordinates {(0,.35)(2,.35)};
\draw[<->] (axis cs:.8,.85)--(axis cs:1.2,.85);
\draw[<->] (axis cs:0,.35)--(axis cs:2,.35);
\draw[dashed] (axis cs:1,.1)--(axis cs:1,1.1);
\node[tinytext,anchor=south] at (axis cs:1,1.08){$\bar u=1$};
\node[tinytext,anchor=west] at (axis cs:1.25,.9){$2\times$ at each point};
\node[tinytext,anchor=west] at (axis cs:.1,.16){$N=4$, identical noise variance};\end{axis}
\node[eq] at (3.9,9.62){$y_k=v+u_kr+w_k,\quad w_k\overset{\mathrm{iid}}\sim\mathcal N(0,\sigma^2)$};
\draw[arr] (9.1,11.0)--(15.9,11.0) node[above left] {$\spanop\{u\}$};
\draw[arr] (9.6,11.0)--(14.0,14.4) node[above] {$\one$};
\draw[dashed,wire] (14,11)--(14,14.4);
\draw[arr,line width=.95pt] (14,11)--(14,14.4) node[right,midway] {$M_u\one$};
\draw[wire] (13.8,11)--(13.8,11.2)--(14,11.2);
\draw (11,11) arc[start angle=0,end angle=37.7,radius=1.4];
\node at (11.45,11.45){$\phi$};
\node[eq] at (12.35,10.05){$M_u=I-uu\trans/(u\trans u),\quad d_0=\sin^2\phi$};
'''
    b+=panel(.1,9.15,7.6,'(c) Exact profiled information')+panel(8.35,9.15,8.05,'(d) Information lost as the angle closes')
    b+=r'''\node[eq] at (3.9,8.32){$\rho=\frac{\one\trans u}{\sqrt{N}\|u\|},\quad d_0=1-\rho^2=1-\bar u^2/\overline{u^2}$};
\node[eq] at (3.9,7.27){$J_v=\sigma^{-2}\one\trans M_u\one=Nd_0/\sigma^2$\\[5pt]
$\widehat v=\frac{\one\trans M_uy}{Nd_0},\quad\Var(\widehat v)=\frac{\sigma^2}{Nd_0}$};
\node[eq] at (3.9,5.71){$\begin{array}{c|c|c|c}\text{Design}&d_0&\rho^2&\Var(\widehat v)/\sigma^2\\\hline
\text{Narrow}&1/26&25/26&13/2\\\text{Wide}&1/2&1/2&1/2\end{array}$};
\node[eq] at (3.9,4.25){$M_u\one=\begin{cases}(3,3,-2,-2)\trans/13&\text{narrow},\\(1,1,0,0)\trans&\text{wide}.\end{cases}$};
\node[tinytext,align=center] at (3.9,2.85){At the same $N$ and $\sigma$, wide excitation gives\\13 times more information about $v$.\\More samples at one fixed nonzero current do not\\separate $v$ from the unknown slope $r$.};
\begin{semilogyaxis}[mathaxis,at={(9.05cm,3.75cm)},anchor=south west,width=7.15cm,height=5.0cm,
 xmin=0,xmax=90,ymin=1,ymax=1000,xtick={0,15,30,45,60,90},ytick={1,10,100,1000},
 xlabel={Column angle $\phi$ [degrees]},ylabel={$N\Var(\widehat v)/\sigma^2$}]
\addplot[domain=2:90,samples=150,line width=.85pt]{1/(sin(x)^2)};
\addplot[only marks,mark=o,mark size=2.2pt] coordinates {(11.309932,26)};
\addplot[only marks,mark=*,mark size=2pt] coordinates {(45,2)};
\node[tinytext,anchor=south west] at (axis cs:13,26){Narrow};
\node[tinytext,anchor=south west] at (axis cs:47,2){Wide};\end{semilogyaxis}
\node[eq] at (12.35,2.45){$u=u_0\one\ (u_0\neq0)\Rightarrow M_u\one=0$\\[5pt]
$\ker[\one\ u]=\spanop\{(-u_0,1)\trans\}$};
\draw[wire] (.1,1.42)--(16.4,1.42);
\node[tinytext,align=center] at (8.25,.78){Design chain: current spread $\longrightarrow$ column angle $\longrightarrow$ retained information $\longrightarrow$ target precision.\\This is a one-path subproblem; multi-path and shared dead-time designs still require the full support test.};
\end{tikzpicture}'''
    add('07_excitation_geometry','Excitation is an angle: visualize information before collecting more samples',b,
        r'''\textbf{Figure 7.} Analytic two-regressor experiment, with four samples and an unknown resistance slope. Duplicate marks denote repeated inputs, not hidden data. The observation-space sketch illustrates an orthogonal decomposition; numerical angles and variance factors are calculated exactly from the displayed designs. The zero-current-only design needs separate handling because the slope column vanishes. No hardware trajectory or operating recommendation is inferred from this normalized example.''')

    b=picture()+panel(.1,15.45,7.6,'(a) A likelihood tube, not a bounded ellipsoid')+panel(8.35,15.45,8.05,'(b) A point estimate along a free fibre')
    # Oblique affine view: a1=(1,0), a2=(.14,.12), b=(.25,.8).
    b+=r'''\begin{scope}[shift={(3.0,11.65)},scale=.75]
\draw[arr] (-1.5,0)--(2.1,0) node[right] {$a_1$};
\draw[arr] (-.85,-.73)--(.9,.77) node[above right] {$a_2$};
\draw[arr] (-.6,-1.92)--(1.3,4.16) node[above] {$b$};
\draw[dashed,wire] plot[domain=0:360,samples=100] ({cos(\x)+.7*sin(\x)},{.6*sin(\x)});
\draw[line width=.9pt] plot[domain=0:360,samples=100] ({.75+cos(\x)+.7*sin(\x)},{2.4+.6*sin(\x)});
\draw[line width=.9pt] plot[domain=0:360,samples=100] ({-.375+cos(\x)+.7*sin(\x)},{-1.2+.6*sin(\x)});
\draw[wire] (-1.596,-1.544)--(-.471,2.056);
\draw[wire] (.846,-.856)--(1.971,2.744);
\draw[densely dotted,->] (1.971,2.744)--(2.284,3.744);
\draw[densely dotted,->] (-1.596,-1.544)--(-1.75,-2.04);
\node[tinytext,align=left,anchor=west] at (2.0,1.5){Unbounded\\null coordinate};
\end{scope}
\node[eq] at (3.9,9.53){$\mathcal T=\{(a_1,a_2,b):a_1^2+a_2^2/25\leq1,\ b\in\R\}$};
\draw[arr] (9.3,12.35)--(15.85,12.35) node[below] {$b$};
\draw[line width=.95pt] (9.65,12.35)--(15.45,12.35);
\fill (11.0,12.35) circle(2pt);\draw[fill=white] (13.95,12.35) circle(2.2pt);
\node[above] at (11,12.52){$\widehat\theta_0$};\node[above] at (13.95,12.52){$\widehat\theta_0+V_0c$};
\draw[<->] (11,11.7)--(13.95,11.7) node[midway,below] {$c$};
\node[eq] at (12.35,14.42){$Z_{\mathrm{eff}}=\diag(1,1/5,0)$\\[5pt]
$\mathcal L(\theta+V_0c;y)=\mathcal L(\theta;y),\quad c\in\R$};
\node[tinytext,align=center] at (12.35,10.36){Pseudoinverse selects $b=0$ on this fibre.\\A prior or a constraint can select another point;\\the measurement likelihood itself stays flat in $b$.};
'''
    b+=panel(.1,8.95,7.6,'(c) Read target variance from tangent lines')+panel(8.35,8.95,8.05,'(d) A tiny null component changes estimability')
    b+=r'''\begin{axis}[mathaxis,axis equal image,at={(1.65cm,4.15cm)},anchor=south west,width=6.5cm,height=5.1cm,
 xmin=-6,xmax=6,ymin=-6,ymax=6,xtick={-5,0,5},ytick={-5,0,5},xlabel={$a_1$},ylabel={$a_2$}]
\addplot[domain=0:360,samples=120,line width=.9pt] ({cos(x)},{5*sin(x)});
\addplot[domain=-6:6,samples=2,dashed]{sqrt(26)-x};
\addplot[domain=-6:6,samples=2,dashed]{-sqrt(26)-x};
\addplot[only marks,mark=*,mark size=1.8pt] coordinates {(.196116,4.902903)(-.196116,-4.902903)};
\draw[->] (axis cs:0,0)--(axis cs:2,2);
\node[tinytext,anchor=west] at (axis cs:2,2){$l$};\end{axis}
\node[eq] at (3.9,2.35){$l=(1,1,0)$\\[4pt]
$\max_{a\trans Fa\leq1}|l\trans a|=\sqrt{l\trans F^+l}=\sqrt{26}$\\[5pt]
$a_*=(1,25,0)\trans/\sqrt{26}$ (with $b=0$)};
\draw[arr] (10,5.95)--(15.35,5.95) node[below] {$a_1$};
\draw[arr] (10,5.95)--(10,8.35) node[left] {$b$};
\draw[arr,line width=.95pt] (10,5.95)--(14.35,5.95) node[below] {$l_0=e_1$};
\draw[arr,dashed] (10,5.95)--(14.35,7.15) node[above] {$l_\psi$};
\draw[wire] (11.35,5.95) arc[start angle=0,end angle=15.4,radius=1.35];
\node[tinytext] at (11.68,6.27){$\psi$};
\node[eq] at (12.35,4.91){$l_\psi=\cos\psi\,e_1+\sin\psi\,e_3$\\[5pt]
$l_\psi\trans(\theta+ce_3)-l_\psi\trans\theta=c\sin\psi$};
\node[eq] at (12.35,3.48){$\sin\psi\neq0\Rightarrow\sup_{\theta\in\mathcal T}|l_\psi\trans\theta|=\infty$\\[5pt]
$l_\psi\trans F^+l_\psi=\cos^2\psi\quad\text{is not a valid bound.}$};
\node[tinytext,align=center] at (12.35,2.03){Support is an exact geometric condition.\\A nearly supported target still needs a justified\\bound or independent measurement of the null part.};
\draw[wire] (.1,1.15)--(16.4,1.15);
\node[tinytext,align=center] at (8.25,.57){Panel (a) uses an oblique view to expose the cylinder; panel (c) uses equal physical axis scales.\\Unit Mahalanobis sets show geometry. A confidence level requires the appropriate quantile and dimension.};
\end{tikzpicture}'''
    add('08_likelihood_fibre','See the quotient: a bounded target inside an unbounded likelihood tube',b,
        r'''\textbf{Figure 8.} Exact normalized singular example from Figure 3. The tube is centered at the fitted identifiable coordinates, translated to the origin for drawing. Tangent lines in (c) display the dual-norm support calculation for an estimable contrast. Any nonzero null coefficient in (d) gives unbounded variation on the same likelihood fibre. The cylinder view is schematic in perspective, with finite displayed length and continuation arrows; it is mathematically unbounded.''')

    b=picture()+panel(.1,15.45,7.6,'(a) Temperature truth and a lagged record')+panel(8.35,15.45,8.05,'(b) The corrected signal can imitate health drift')
    b+=r'''\begin{axis}[mathaxis,at={(.85cm,10.55cm)},anchor=south west,width=7.15cm,height=4.55cm,
 xmin=-.5,xmax=4,ymin=-.5,ymax=9,xtick={0,1,2,3,4},ytick={0,4,8},xlabel={$t/\tau$},ylabel={Temperature rise [K]},
 legend style={at={(.97,.06)},anchor=south east}]
\addplot[line width=.85pt] coordinates {(-.5,0)(0,0)(0,8)(4,8)};\addlegendentry{Junction truth}
\addplot[dashed,domain=0:4,samples=100,line width=.85pt]{8*(1-exp(-x))};\addlegendentry{Recorded indication}
\draw[<->] (axis cs:1,5.056964)--(axis cs:1,8) node[midway,left,tinytext] {$-e(t)$};\end{axis}
\node[eq] at (3.9,9.45){$\tau\dot T^{\mathrm{record}}+T^{\mathrm{record}}=T^{\mathrm{true}}$};
\begin{axis}[mathaxis,at={(9.05cm,10.55cm)},anchor=south west,width=7.15cm,height=4.55cm,
 xmin=0,xmax=4,ymin=0,ymax=27,xtick={0,1,2,3,4},ytick={0,5,10,20,25},xlabel={$t/\tau$},ylabel={Corrected change [mV]},
 legend style={at={(.97,.97)},anchor=north east}]
\addplot[domain=0:4,samples=100,line width=.85pt]{20*exp(-x)};\addlegendentry{Zero health change}
\addplot[domain=0:4,samples=100,dashed,line width=.85pt]{5+20*exp(-x)};\addlegendentry{5 mV health change}
\addplot[domain=0:4,samples=2,densely dotted,forget plot]{5};
\draw[densely dotted] (axis cs:1.386294,0)--(axis cs:1.386294,5);
\addplot[only marks,mark=o,mark size=2.2pt,forget plot] coordinates {(1.386294,5)};\end{axis}
\node[tinytext,align=center] at (12.35,9.45){Analytic step: $\Delta T=8$ K, $\alpha=2.5$ mV/K.};
'''
    b+=panel(.1,9.0,7.6,'(c) Derive the sign and the settling budget')+panel(8.35,9.0,8.05,'(d) Visual error budget for an estimable target')
    b+=r'''\node[eq] at (3.9,8.1){$e(t)=T^{\mathrm{record}}-T^{\mathrm{true}}=-\Delta T e^{-t/\tau}$};
\node[eq] at (3.9,6.97){$\widehat{\Delta h}(t)=\Delta h-\alpha[e(t)-e_b]$\\[5pt]
$e_b=0\Rightarrow\widehat{\Delta h}(t)=\Delta h+20e^{-t/\tau}\ \mathrm{mV}$};
\node[eq] at (3.9,5.42){$|\alpha\Delta T|e^{-t/\tau}\leq B_T$\\[5pt]
$t\geq\tau\max\{0,\log(|\alpha\Delta T|/B_T)\}$};
\node[eq] at (3.9,4.26){$B_T=5\ \mathrm{mV}\Rightarrow t\geq\tau\log4\simeq1.386\tau$};
\node[tinytext,align=center] at (3.9,2.86){The threshold controls this thermal residual only.\\It does not by itself control random false alarms.\\The step and first-order lag are illustrative;\\a real lag bound needs sensor or device evidence.};
\draw[arr] (8.8,5.55)--(16.0,5.55) node[below] {$\Delta h$};
\draw[densely dotted,wire] (12.4,5.55)--(12.4,6.0);
\draw[line width=1.2pt] (10.6,6.0)--(14.2,6.0);
\draw[wire] (10.6,5.83)--(10.6,6.17);\draw[wire] (14.2,5.83)--(14.2,6.17);
\fill (12.4,6) circle(2pt);\node[above] at (12.4,6.18){$\widehat{\Delta h}$};
\draw[<->] (10.6,6.85)--(14.2,6.85) node[midway,above] {$2zs_\Delta$};
\draw[dashed,line width=.8pt] (9.6,6)--(10.6,6);\draw[dashed,line width=.8pt] (14.2,6)--(15.2,6);
\draw[wire] (9.6,5.7)--(9.6,6.3);\draw[wire] (15.2,5.7)--(15.2,6.3);
\draw[<->] (14.2,5.2)--(15.2,5.2) node[midway,below] {$B_T$};
\node[eq] at (12.35,8.06){$\widehat{\Delta h}=\Delta h+b_T+s_\Delta Z,\quad Z\sim\mathcal N(0,1)$};
\node[eq] at (12.35,4.01){$|b_T|\leq B_T,\quad z=\Phi^{-1}(1-\alpha_0/2)$\\[5pt]
$\mathcal I=[\widehat{\Delta h}-(zs_\Delta+B_T),\ \widehat{\Delta h}+(zs_\Delta+B_T)]$};
\node[tinytext,align=center] at (12.35,2.72){For known $s_\Delta$ and a valid deterministic bias bound,\\$\Pr(\Delta h\in\mathcal I)\geq1-\alpha_0$.\\Solid: Gaussian radius; dashed: bounded-bias extension.};
\draw[wire] (.1,1.35)--(16.4,1.35);
\node[tinytext,align=center] at (8.25,.72){Temperature correction $\longrightarrow$ residual lag $\longrightarrow$ target bias $\longrightarrow$ uncertainty allocation.\\Averaging windows and time-varying temperature require the complete estimator's weighted discrepancy operator.};
\end{tikzpicture}'''
    add('09_thermal_dynamics','Thermal lag: a transient bias can look like a persistent parameter change',b,
        r'''\textbf{Figure 9.} Pointwise analytic illustration, without noisy samples or simulated sensor evidence. The baseline is settled and the junction temperature undergoes a step; the indication follows a first-order lag. Curves compare zero health change and a true 5 mV change. The interval in (d) adds a deterministic bias budget to Gaussian uncertainty; it is valid only for an estimable target, correct variance, and a justified bound. It is not an inferred uncertainty statement for the retained SIL trials.''')

    b=picture()+panel(.1,15.45,7.6,'(a) One reference, correlated differences')+panel(8.35,15.45,8.05,'(b) A per-test rate accumulates under repetition')
    b+=r'''\node[device,minimum width=1.75cm] (ref) at (1.6,13.18){Reference $R$};
\node[device,minimum width=1.55cm] (m1) at (5.65,14.35){Record 1};
\node[device,minimum width=1.55cm] (m2) at (5.65,13.18){Record 2};
\node[device,minimum width=1.55cm] (mm) at (5.65,12.01){Record $m$};
\draw[arr] (ref.east)--(m1.west);\draw[arr] (ref.east)--(m2.west);\draw[arr] (ref.east)--(mm.west);
\node[tinytext] at (5.65,12.6){$\vdots$};
\node[eq] at (3.9,10.69){$\epsilon_i=\xi_i-\xi_R,\quad\Cov(\epsilon)=vI_m+v_R\one\one\trans$\\[5pt]
$v_R=v\Rightarrow\operatorname{Corr}(\epsilon_i,\epsilon_j)=1/2\quad(i\neq j)$};
\begin{axis}[mathaxis,at={(9.05cm,10.5cm)},anchor=south west,width=7.15cm,height=4.6cm,
 xmin=1,xmax=60,ymin=0,ymax=1.05,xtick={1,10,20,40,60},ytick={0,.5,1},xlabel={Number of tests $m$},ylabel={At least one null alert},
 legend style={at={(.97,.05)},anchor=south east}]
\addplot[domain=1:60,samples=100,line width=.85pt]{1-.95^x};\addlegendentry{Independent tests}
\addplot[domain=1:20,samples=40,dashed,line width=.85pt]{.05*x};
\addplot[domain=20:60,samples=2,dashed,line width=.85pt,forget plot]{1};\addlegendentry{Union bound}
\addplot[only marks,mark=o,mark size=2pt,forget plot] coordinates {(20,.641514)};
\draw[wire] (axis cs:20,.641514)--(axis cs:23,.58);
\node[tinytext,anchor=north west] at (axis cs:23,.57){$64.2\%$ at $m=20$};\end{axis}
\node[tinytext,align=center] at (12.35,9.45){Independence is an explicit scenario, not implied by (a).};
'''
    b+=panel(.1,9.05,7.6,'(c) State exactly which error rate is controlled')+panel(8.35,9.05,8.05,'(d) Calibrate the whole task if that is the claim')
    b+=r'''\node[eq] at (3.9,8.08){$\Pr(A_i)\leq\alpha_0\ \ \forall i$\\[5pt]
$\Pr\!\left(\bigcup_{i=1}^m A_i\right)\leq\min\{1,m\alpha_0\}$};
\node[eq] at (3.9,6.7){$\text{Independent: }\Pr\!\left(\bigcup_i A_i\right)=1-(1-\alpha_0)^m$};
\node[eq] at (3.9,5.4){$n_{\mathrm{cal}}=200,\ k=191:\quad\alpha_{\mathrm{rank}}=10/201$\\[5pt]
$\Pr_{\mathrm{cal,new}}(S_{\mathrm{new}}>S_{(191)})\leq\alpha_{\mathrm{rank}}$};
\node[tinytext,align=center] at (3.9,3.9){This is marginal over the calibration split and new case.\\The same threshold used across many tests can\\couple decisions, even with independent raw records.};
\node[tinytext,align=center] at (3.9,2.59){Per-target, per-ADC calibration in Figure 6 does not\\establish a trajectory-level or multi-target guarantee.\\A guarantee must match the unit of exchangeability.};
\node[draw,align=center,minimum width=6.8cm] (traj) at (12.35,8.05){Freeze horizon $m$, targets $j=1,\ldots,J$,\\estimator, reference rule and variance model};
\node[draw,align=center,minimum width=6.8cm] (score) at (12.35,6.42){One score per complete null trajectory\\$S_i^{\max}=\max_{t\leq m,\,j\leq J}|\widehat{\Delta h}_{i,t,j}|/s_{i,t,j}$};
\node[draw,align=center,minimum width=6.8cm] (quant) at (12.35,4.78){$q_{\max}=S_{(191)}^{\max}$ from 200 trajectories\\Alert if $S_{\mathrm{new}}^{\max}>q_{\max}$};
\draw[arr] (traj.south)--(score.north);\draw[arr] (score.south)--(quant.north);
\node[eq] at (12.35,3.13){$\Pr_{\mathrm{cal,new}}\{\text{any null alert in the frozen task}\}$\\[5pt]
$\leq10/201$ under trajectory exchangeability.};
\node[tinytext,align=center] at (12.35,1.95){Proposed validation protocol; not evaluated in existing trials.};
\draw[wire] (.1,1.3)--(16.4,1.3);
\node[tinytext,align=center] at (8.25,.66){Measurement claim $\longrightarrow$ dependence model $\longrightarrow$ calibration unit $\longrightarrow$ reported error rate.\\Exchangeability of complete trajectories permits within-trajectory dependence; adaptive horizons need new justification.};
\end{tikzpicture}'''
    add('10_repeated_decisions','From one calibrated decision to a complete monitoring task',b,
        r'''\textbf{Figure 10.} Reference reuse gives an explicit dependence mechanism. Panel (b) shows the exact independent-test familywise rate at $\alpha_0=0.05$ alongside the dependence-agnostic union bound, rather than observed alarm data. Panel (d) is a proposed finite-horizon calibration protocol: exchangeable complete null trajectories are its sampling units. Its rank guarantee is marginal and requires the entire task to be frozen before calibration. Existing per-target SIL rates are retained as narrower evidence.''')
