"""Monochrome mathematical plates; exact fixtures and analytic plots, no new trials."""
from pathlib import Path
import json, hashlib, re, math
import numpy as np
from scipy.stats import norm
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'figures/latex_illustrations'; DEST.mkdir(parents=True,exist_ok=True)
DATA=json.loads((ROOT/'data/independent_sil_results.json').read_text())
PRE=r'''\documentclass[10pt]{article}
\usepackage[T1]{fontenc}
\usepackage{lmodern,amsmath,amssymb,booktabs,array,graphicx,tikz,pgfplots}
\usepackage[paperwidth=190mm,paperheight=224mm,margin=10mm]{geometry}
\usepackage[hidelinks]{hyperref}
\hypersetup{pdfauthor={},pdftitle={Monochrome mathematical illustrations: semiconductor-health identifiability}}
\usetikzlibrary{arrows.meta,positioning,calc,patterns}
\pgfplotsset{compat=1.18}
\newcommand{\trans}{^{\mathsf T}}\newcommand{\R}{\mathbb R}
\newcommand{\one}{\mathbf 1}\newcommand{\zero}{\mathbf 0}
\DeclareMathOperator{\rank}{rank}\DeclareMathOperator{\range}{range}
\DeclareMathOperator{\diag}{diag}\DeclareMathOperator{\Cov}{Cov}
\DeclareMathOperator{\Var}{Var}\DeclareMathOperator{\spanop}{span}
\newcommand{\FigureNote}{Visualization code prepared with assistance from OpenAI Codex.}
\pagestyle{empty}\setlength{\parindent}{0pt}
\tikzset{
 every picture/.style={font=\fontsize{9}{11}\selectfont,>=Latex},
 every node/.style={font=\fontsize{9}{11}\selectfont,inner sep=2pt},
 head/.style={font=\fontsize{9.5}{12}\selectfont\bfseries,anchor=west},
 tinytext/.style={font=\fontsize{8}{10}\selectfont},eq/.style={align=center},
 wire/.style={draw=black,line width=.5pt},arr/.style={->,black,line width=.65pt},
 device/.style={draw,fill=white,minimum width=1.1cm,minimum height=.38cm}}
\pgfplotsset{mathaxis/.style={axis lines=left,axis line style={black,line width=.5pt},
 tick style={black},tick label style={font=\fontsize{8}{10}\selectfont},
 label style={font=\fontsize{8.5}{10}\selectfont},grid=major,
 major grid style={black!15,line width=.2pt},
 legend style={font=\fontsize{8}{10}\selectfont,draw=none,fill=white},
 every axis plot/.append style={black}}}
\newcommand{\PlateTitle}[2]{\small\textbf{Mathematical plate #1}\hfill Monochrome TikZ / PGFPlots\par
 \vspace{1.5mm}{\large\bfseries #2}\par\vspace{2mm}}
\newcommand{\PlateCaption}[1]{\par\vspace{1mm}{\fontsize{8.5}{10.5}\selectfont #1\par
 \vspace{1mm}\FigureNote\par}}
'''
figures=[]
def add(slug,title,body,caption):figures.append(dict(slug=slug,title=title,body=body,caption=caption))
def picture():return r'\begin{tikzpicture}[x=1cm,y=1cm]\path[use as bounding box] (0,0) rectangle (16.6,15.8);'+'\n'
def panel(x,y,w,label):
    return r'\node[head] at ('+f'{x},{y}'+r'){'+label+r'};\draw[wire] ('+f'{x},{y-.28}'+')--('+f'{x+w},{y-.28}'+');\n'

b=picture()+panel(.1,15.45,7.6,'(a) Paths define the linear observation')+panel(8.35,15.45,8.05,'(b) Exact incidence in an explicit ordering')
paths=['P^+','N^+','P^-','N^-','O^+','O^-']
devices=['S_1','D_4','D_1','S_4','S_2','D_3','S_3','D_2']
yp=[14.4-.72*j for j in range(6)];yd=[14.4-.58*j for j in range(8)]
for i,name in enumerate(paths):b+=r'\node[anchor=east] (p'+str(i)+') at (1.9,'+str(yp[i])+r'){$'+name+r'$};'+'\n'
for j,name in enumerate(devices):b+=r'\node[device] (d'+str(j)+') at (5.7,'+str(yd[j])+r'){$'+name+r'$};'+'\n'
for i,j in [(0,0),(1,1),(2,2),(3,3),(4,4),(4,5),(5,6),(5,7)]:b+=r'\draw[arr] (p'+str(i)+'.east)--(d'+str(j)+'.west);\n'
b+=r'''\node[tinytext,align=left,anchor=west] at (.5,9.7){$P,O,N$: pole levels; superscript: current sign.\\
Singletons reveal outer-device coefficients;\\neutral paths reveal two series sums.};
\node[eq] at (12.35,14.65){$\vartheta=(S_1,D_4,D_1,S_4,S_2,D_3,S_3,D_2)\trans$};
\node[eq] at (12.35,12.85){$A=\left[\begin{array}{cccc|cccc}
1&0&0&0&0&0&0&0\\0&1&0&0&0&0&0&0\\0&0&1&0&0&0&0&0\\0&0&0&1&0&0&0&0\\
0&0&0&0&1&1&0&0\\0&0&0&0&0&0&1&1
\end{array}\right],\quad p=A\vartheta$};
\node[eq] at (12.35,10.95){$\zeta=(1,1,-1,-1,1,-1)\trans$\\[5pt]
$\rank A=6,\quad\dim\ker A=8-6=2$};
'''
b+=panel(.1,8.65,7.6,'(c) Kernel and canonical representative')+panel(8.35,8.65,8.05,'(d) Two free kernel coordinates')
b+=r'''\node[eq] at (2.25,6.65){$K=\begin{bmatrix}
0&0\\0&0\\0&0\\0&0\\1&0\\-1&0\\0&1\\0&-1\end{bmatrix}$};
\node[eq] at (5.65,6.85){$AK=0,\quad K\trans K=2I_2$\\[8pt]
$P_A=A^+A$\\[4pt]$=\diag(I_4,\tfrac12J_2,\tfrac12J_2)$\\[5pt]
$J_2=\begin{bmatrix}1&1\\1&1\end{bmatrix}$};
\node[eq] at (3.9,4.9){$\vartheta^\dagger=A^+p,\qquad\vartheta=\vartheta^\dagger+K\binom{s}{t}$};
\node[eq] at (3.9,3.72){$\begin{aligned}
\vartheta_5&=p_5/2+s,&\vartheta_6&=p_5/2-s,\\
\vartheta_7&=p_6/2+t,&\vartheta_8&=p_6/2-t.\end{aligned}$};
\node[tinytext,align=center] at (3.9,2.1){The minimum-norm split is a coordinate convention.\\It supplies no new device-level measurement.};
\draw[arr] (9.0,3.7)--(15.9,3.7) node[below] {$s$};
\draw[arr] (12.0,1.9)--(12.0,7.65) node[left] {$t$};
\draw[dashed,wire] (9.55,2.4) rectangle (15.1,7.0);
\fill (12,3.7) circle(1.8pt);\node[anchor=north east] at (11.9,3.55){$\vartheta^\dagger$};
\fill (13.9,5.5) circle(1.8pt);\node[anchor=west] at (14.0,5.5){$\vartheta_1$};
\draw[fill=white] (10.3,6.3) circle(2pt);\node[anchor=east] at (10.2,6.3){$\vartheta_2$};
\draw[<->,wire] (10.3,6.3)--(13.9,5.5);
\node[tinytext,align=center] at (12.45,7.95){Every point gives the same $A\vartheta$.};
\node[tinytext,align=center] at (12.45,1.5){Kernel coordinates $(s,t)$; the dashed frame\\is a drawing boundary, not a physical constraint.};
\draw[wire] (.1,.96)--(16.4,.96);
\node[eq] at (8.25,.5){$l\trans\vartheta\text{ is determined by }p\ \Longleftrightarrow\ K\trans l=0
\ \Longleftrightarrow\ l_5=l_6,\ l_7=l_8.$};\end{tikzpicture}'''
add('01_incidence_kernel','Conduction incidence, null directions and the observable quotient',b,
r'''\textbf{Figure 1.} Exact T-type incidence analysis for one health channel in the ideal per-leg pole-drop observation. Device symbols in $\vartheta$ denote coefficient coordinates. The reordered columns expose the two inseparable neutral pairs. Panels (c)--(d) show the entire unrestricted inverse image, rather than a single fitted solution. The same incidence ambiguity applies separately to threshold and resistance channels; the complete three-phase current model imposes additional restrictions.''')

b=picture()+panel(.1,15.45,7.6,'(a) Complete observable current operator')+panel(8.35,15.45,8.05,'(b) Isolated-neutral common-mode projection')
b+=r'''\node[eq] at (3.9,14.36){$C=\begin{bmatrix}1/\sqrt2&-1/\sqrt2&0\\1/\sqrt6&1/\sqrt6&-2/\sqrt6\end{bmatrix}$};
\node[eq] at (3.9,13.32){$CC\trans=I_2,\quad C\one=0,\quad C\trans C=I_3-\tfrac13\one\one\trans$};
\node[eq] at (3.9,12.25){$\begin{aligned}B_k&=\diag(a_{p_a(k)}\trans,a_{p_b(k)}\trans,a_{p_c(k)}\trans),\\
Z_k&=[-CD_{z,k}B_k\quad-CD_{i,k}B_k],\\H_k&=[-C(n_k\odot z_k)\quad-Ci_k].\end{aligned}$};
\node[eq] at (3.9,10.85){$X=\operatorname{stack}_k[Z_k\ H_k]\in\R^{2N\times50}$\\[4pt]
$x=(v_a,v_b,v_c,r_a,r_b,r_c,\beta_{\mathrm{dt}},\Delta R_s)$};
\draw[wire] (9.0,11.3)--(12.0,10.4)--(15.6,11.45)--(12.6,12.35)--cycle;
\node[tinytext] at (14.4,10.8){$\one^\perp$};
\draw[arr] (12.3,11.25)--(12.3,14.5) node[above] {$\spanop\{\one\}$};
\draw[dashed,wire] (11.2,11.65)--(11.2,14.1);
\fill (11.2,12.6) circle(2pt);\node[anchor=east] at (11.1,12.6){$\delta u$};
\draw[fill=white] (11.2,14.1) circle(2pt);\node[anchor=west] at (11.35,14.1){$\delta u+s\one$};
\fill (11.2,11.65) circle(2pt);\node[anchor=east,tinytext] at (11.05,11.6){$P_0\delta u$};
\node[eq] at (12.3,10.05){$P_0=C\trans C,\qquad C(\delta u+s\one)=C\delta u$};
'''
b+=panel(.1,9.3,16.3,'(c) Exact kernel under finite separating conditions')
b+=r'''\node[eq] at (8.25,8.55){$\ker X=\underbrace{(\ker A)^6\times\{(0,0)\}}_{\mathcal K_A:\;12\text{ dimensions}}
\ \oplus\ \underbrace{\spanop\{t_v,t_r\}}_{2\text{ common gauges}}$};
\node[eq] at (4.1,7.55){$h=(-1,-1,1,1,-1,0,1,0)\trans$\\[4pt]$Ah=-\zeta$};
\node[eq] at (12.3,7.55){$g=(1,1,1,1,1,0,1,0)\trans$\\[4pt]$Ag=\one$};
\node[eq] at (8.25,6.55){$\begin{aligned}t_v&=(h,h,h,0,0,0,0,0),\\t_r&=(0,0,0,g,g,g,0,-1).
\end{aligned}\qquad\rank X=50-(12+2)=36$};
\node[tinytext,align=center] at (8.25,5.46){All paths visited; connected co-occurrence graph; affine spanning of balanced currents;\\a repeated path/current observation with $C[(n-n')\odot z]\neq0$. These are sufficient conditions.};
\node[eq] at (8.25,4.7){$n_x(k)\equiv\bar n\neq0:\quad t_{\mathrm{dt}}=(-\bar n g,-\bar n g,-\bar n g,0,0,0,1,0),\quad\rank X=35.$};
'''
b+=panel(.1,3.92,16.3,'(d) Annihilate the entire kernel before naming a measurand')
b+=r'''\node[eq] at (4.1,2.93){$\begin{aligned}
h_v&=v_{S1a}-v_{D4a}:&(e_1-e_2)\trans h&=0,\\
h_r&=r_{S1a}-r_{S4a}:&(e_1-e_4)\trans g&=0.\end{aligned}$};
\node[eq] at (12.3,2.93){$v_{S1a}:\ e_1\trans h=-1\neq0$\\[5pt]Absolute threshold attribution fails.};
\node[tinytext,align=center] at (8.25,1.72){The contrasts also annihilate the incidence directions in their health channel.\\Support here assumes the specified nuisance model and separating excitation.};
\draw[wire] (.1,1.08)--(16.4,1.08);
\node[eq] at (8.25,.55){Unknown scale adds $h_v'=\lambda h_v+(1-\lambda)(V_{\mathrm{dc}}-h_v^0)$; an independent machine-scale reference is still needed.};
\end{tikzpicture}'''
add('02_joint_kernel','From path ambiguity to the complete three-phase kernel',b,
r'''\textbf{Figure 2.} Common-mode and stator-resistance gauges are additional to the six incidence kernels. The vectors $h,g$ use Figure 1's ordering; substitution verifies both gauge directions exactly. Projection geometry is schematic, while the matrices and rank counts are exact. Rank 36 is conditional on the theorem's separating conditions, not guaranteed for every controller trajectory. Fixed commutation adds another null direction; unknown machine scale introduces a further model-specific equivalence.''')

b=picture()+panel(.1,15.45,7.6,'(a) Whiten and remove nuisance')+panel(8.35,15.45,8.05,'(b) Identifiable and null SVD coordinates')
b+=r'''\draw[arr] (.9,11.3)--(7.15,11.3) node[above left] {$\range(WH)$};
\draw[arr] (1.15,11.3)--(5.05,14.1) node[above] {$WZt$};
\draw[dashed,wire] (5.05,11.3)--(5.05,14.1);
\draw[arr] (1.15,11.3)--(5.05,11.3);
\draw[arr,line width=.9pt] (5.05,11.3)--(5.05,14.1) node[midway,right] {$MWZt$};
\draw[wire] (4.85,11.3)--(4.85,11.5)--(5.05,11.5);
\node[eq] at (3.9,10.2){$\Pi_{WH}=WH(WH)^+,\quad M=I-\Pi_{WH}$\\[4pt]$M=M\trans=M^2,\quad MWH=0$};
\node[eq] at (12.35,14.45){$Z_{\mathrm{eff}}=MWZ=U_rD_rV_r\trans$\\[5pt]$D_r=\diag(\sigma_1,\ldots,\sigma_r),\quad\sigma_j>0$};
\node[eq] at (12.35,13.1){$\theta=V_ra+V_0b,\quad Z_{\mathrm{eff}}V_0=0$\\[5pt]
$\bar y=U_r\trans MW\widetilde y=D_ra+\varepsilon_r$\\[4pt]$\varepsilon_r\sim\mathcal N(0,I_r)$};
\node[eq] at (12.35,11.36){$F_\theta=V_rD_r^2V_r\trans$\\[5pt]$F_\theta^+=V_rD_r^{-2}V_r\trans$};
\node[tinytext,align=center] at (12.35,10.15){White noise has identity covariance on the retained $U_r$\\coordinates, not on the full projected ambient space.};
'''
b+=panel(.1,9.35,7.6,'(c) Estimability precedes precision')+panel(8.35,9.35,8.05,'(d) Small singular value, large variance')
b+=r'''\node[eq] at (3.9,8.3){$l=V_r\ell_r+V_0\ell_0$\\[5pt]$l\in\range(F_\theta)\ \Longleftrightarrow\ \ell_0=0$};
\node[eq] at (3.9,6.86){$\widehat h=\ell_r\trans D_r^{-1}\bar y$\\[6pt]$\Var(\widehat h)=\sum_{j=1}^{r}\ell_{r,j}^2/\sigma_j^2$};
\node[eq] at (3.9,5.36){$Z_{\mathrm{eff}}=\diag(1,\epsilon,0),\quad\epsilon=0.2$\\[5pt]$F_\theta=\diag(1,0.04,0)$};
\node[eq] at (3.9,3.61){$\begin{array}{c|c|c}l&\text{Estimable?}&\text{Variance}\\\hline
e_1&\text{yes}&1\\e_2&\text{yes}&25\\e_3&\text{no}&\text{undefined}\end{array}$};
\node[tinytext,align=center] at (3.9,2.27){Yet $e_3\trans F_\theta^+e_3=0$.\\That finite number is not an uncertainty bound.};
\begin{axis}[mathaxis,axis equal image,at={(9.1cm,3.1cm)},anchor=south west,width=7.1cm,height=6.5cm,
 xmin=-6,xmax=6,ymin=-6,ymax=6,xtick={-5,0,5},ytick={-5,0,5},
 xlabel={Identifiable error $a_1$},ylabel={Identifiable error $a_2$},clip=false]
\addplot[domain=0:360,samples=121,line width=.85pt] ({cos(x)},{5*sin(x)});
\addplot[only marks,mark=*,mark size=1.6pt] coordinates {(1,0)(0,5)};
\node[anchor=west,tinytext] at (axis cs:.4,4.6){$s_{e_2}=5$};
\node[anchor=south west,tinytext] at (axis cs:1.1,.5){$s_{e_1}=1$};
\node[anchor=west,tinytext] at (axis cs:1.5,-3.4){$a_1^2+a_2^2/25=1$};\end{axis}
\node[tinytext,align=center] at (12.35,1.85){Unit Mahalanobis contour; equal axis scales.\\The null coordinate $b$ is absent from the observable ellipse.};
\draw[wire] (.1,1.15)--(16.4,1.15);
\node[eq] at (8.25,.57){$\operatorname{bias}(\widehat h)=\ell_r\trans D_r^{-1}U_r\trans Mb_w,\quad
\|Mb_w\|_2\leq E\ \Rightarrow\ |\operatorname{bias}(\widehat h)|\leq E\sqrt{l\trans F_\theta^+l}.$};
\end{tikzpicture}'''
add('03_projection_svd','Nuisance projection, singular directions and target-specific uncertainty',b,
r'''\textbf{Figure 3.} Derivation of an estimable Gaussian target through whitening and orthogonal nuisance removal. Here $W\trans W=\Sigma^{-1}$ and $b_w=Wb_{\mathrm{mod}}$. The ellipse displays exact covariance in a normalized example, not physical-drive calibration. Small positive singular values amplify variance; a null target cannot be recovered by a pseudoinverse convention. The bias relation requires estimable $l$ and a justified whole-window discrepancy bound $E$.''')

b=picture()+panel(.1,15.45,7.6,'(a) Rank-one shared-nuisance penalty')+panel(8.35,15.45,8.05,'(b) Common and contrast eigenmodes')
b+=r'''\node[eq] at (3.9,14.25){$y_k=v_{p(k)}+u_kr_{p(k)}+n_k\beta+w_k$\\[5pt]$Q=\diag(q),\quad d_0=1-m_1^2/m_2>0$};
\node[eq] at (3.9,12.65){$\frac{\sigma_w^2}{N}F_v=d_0Q-\frac{\mu^2d_0^2}{s_n^2+\mu^2d_0}qq\trans$\\[8pt]
$\frac{N}{\sigma_w^2}F_v^{-1}=d_0^{-1}Q^{-1}+\frac{\mu^2}{s_n^2}\one\one\trans$};
\node[eq] at (3.9,10.85){$q_1=q_2=\tfrac12,\ d_0=\tfrac14,\ \kappa=\mu^2/s_n^2$\\[5pt]
$\mathcal C=\frac{N}{\sigma_w^2}F_v^{-1}=8I_2+\kappa\one\one\trans$};
\begin{axis}[mathaxis,axis equal image,at={(9.0cm,10.5cm)},anchor=south west,width=7.2cm,height=4.2cm,
 xmin=-10,xmax=10,ymin=-4,ymax=4,xtick={-8,0,8},ytick={-3,0,3},
 xlabel={Common error $a_+$},ylabel={Contrast error $a_-$},legend style={at={(.5,1.03)},anchor=south,legend columns=3}]
\addplot[domain=0:360,samples=121,line width=.8pt] ({sqrt(8)*cos(x)},{sqrt(8)*sin(x)});\addlegendentry{$\kappa=0$}
\addplot[domain=0:360,samples=121,line width=.8pt,dashed] ({4*cos(x)},{sqrt(8)*sin(x)});\addlegendentry{$\kappa=4$}
\addplot[domain=0:360,samples=121,line width=.8pt,densely dotted] ({sqrt(80)*cos(x)},{sqrt(8)*sin(x)});\addlegendentry{$\kappa=36$}\end{axis}
'''
b+=panel(.1,9.2,7.6,'(c) Level and contrast variance limits')+panel(8.35,9.2,8.05,'(d) Eigenvalues and a rank counterexample')
b+=r'''\begin{semilogyaxis}[mathaxis,at={(.85cm,3.8cm)},anchor=south west,width=7.15cm,height=5.0cm,
 xmin=0,xmax=1,ymin=5,ymax=1000,xtick={0,.25,.5,.75,1},ytick={10,100,1000},
 xlabel={$\gamma=\mu^2/(\mu^2+s_n^2)$},ylabel={Normalized target variance},legend style={at={(.04,.97)},anchor=north west}]
\addplot[domain=0:.998,samples=151,line width=.8pt]{8+x/(1-x)};\addlegendentry{Level $v_1$}
\addplot[domain=0:.998,samples=2,dashed,line width=.8pt]{16};\addlegendentry{Difference $v_1-v_2$}
\addplot[domain=0:.998,samples=151,densely dotted,line width=.8pt]{8+2*x/(1-x)};\addlegendentry{Unit common target $e_+\trans v$}\end{semilogyaxis}
\node[tinytext,align=center] at (3.9,2.45){Variance is multiplied by $N/\sigma_w^2$.\\Occupancy and conditional moments remain fixed.};
\node[eq] at (12.35,8.24){$e_+=\tfrac1{\sqrt2}(1,1)\trans,\quad e_-=\tfrac1{\sqrt2}(1,-1)\trans$\\[6pt]
$\mathcal C e_+=(8+2\kappa)e_+,\quad\mathcal C e_-=8e_-$};
\node[eq] at (12.35,6.9){$s_n^2\to0:\quad F_v\to\frac{Nd_0}{\sigma_w^2}(Q-qq\trans)$\\[5pt]$\ker F_v=\spanop\{\one\}$};
\node[eq] at (12.35,5.68){$l\trans\one=0:\quad l\trans F_v^+l=\frac{\sigma_w^2}{Nd_0}\sum_p l_p^2/q_p$};
\node[tinytext,align=center] at (12.35,4.63){The level becomes unidentifiable;\\zero-sum contrasts retain finite bounds.};
\node[eq] at (12.35,3.05){$n=u\in\{1,2\}:\quad\begin{bmatrix}1&1&1\\1&2&2\end{bmatrix}\begin{bmatrix}0\\1\\-1\end{bmatrix}=0$};
\node[tinytext,align=center] at (12.35,1.95){Nevertheless $\rho=3/\sqrt{10}<1$ and $\gamma=0.9<1$.\\Scalar excitation indices alone do not certify full rank.};
\draw[wire] (.1,1.12)--(16.4,1.12);
\node[tinytext,align=center] at (8.25,.5){The closed form needs common conditional moments and $\mathbb E(un)=\mathbb E(u)\mathbb E(n)$.\\For correlated or unequal path statistics, evaluate the complete profiled operator.};\end{tikzpicture}'''
add('04_shared_information','Shared dead time: the covariance eigenmode that actually diverges',b,
r'''\textbf{Figure 4.} Exact balanced two-path calculation with $d_0=1/4$. Unit Mahalanobis contours in (b) use the eigenbasis and units $\sigma_w/\sqrt N$. The rank-one term lengthens only the common coordinate. Panel (c) distinguishes an unnormalized difference from the unit common target. The counterexample violates moment factorization: resistance and dead-time columns coincide despite nondegenerate scalar moment indices. These are statistical design calculations, not measured drive operating curves.''')

b=picture()+panel(.1,15.45,7.6,'(a) Temperature offset is an augmented kernel')+panel(8.35,15.45,8.05,'(b) Scalar cross-section of the thermal fibre')
b+=r'''\node[eq] at (3.9,14.22){$y_{\mathrm{corr}}=Z\theta+Z\Gamma b+H\eta+w$\\[5pt]
$b=T^{\mathrm{true}}-T^{\mathrm{record}},\quad\Gamma=\begin{bmatrix}\diag(\alpha_v)\\\diag(\alpha_r)\end{bmatrix}$};
\node[eq] at (3.9,12.4){$[Z\ Z\Gamma]\binom{-\Gamma c}{c}=0$\\[7pt]$(\theta,b)\mapsto(\theta-\Gamma c,b+c)$};
\node[eq] at (3.9,10.85){$[l\trans\ 0]\binom{-\Gamma c}{c}=-l\trans\Gamma c$\\[5pt]Unrestricted $b$ requires $l\trans\Gamma=0$.};
\begin{axis}[mathaxis,at={(9.0cm,10.6cm)},anchor=south west,width=7.2cm,height=4.1cm,
 xmin=0,xmax=20,ymin=-6,ymax=6,xtick={0,10,20},ytick={-5,0,5},
 xlabel={Reference-temperature drift $a$ [mV]},ylabel={Offset $b$ [K]},clip=false]
\addplot[domain=0:20,samples=2,line width=.85pt]{(10-x)/2};
\addplot[domain=0:20,samples=2,dashed,line width=.5pt]{(11-x)/2};
\addplot[domain=0:20,samples=2,dashed,line width=.5pt]{(9-x)/2};
\addplot[only marks,mark=o,mark size=2.2pt] coordinates {(0,5)(10,0)(20,-5)};
\node[tinytext,anchor=west] at (axis cs:1,-4.7){$a+2b=10$ mV};
\node[tinytext,anchor=east] at (axis cs:19,4.7){Same corrected observation};\end{axis}
'''
b+=panel(.1,9.25,7.6,'(c) Linear image of an offset-error box')+panel(8.35,9.25,8.05,'(d) Calibration drift between records')
b+=r'''\node[eq] at (3.9,8.25){$\mathcal E_T=\{-L_{\mathrm{tar}}\Gamma\delta e:\ |\delta e_j|\leq\epsilon_j\}$\\[5pt]
$B_{T,k}=\sum_j |(L_{\mathrm{tar}}\Gamma)_{kj}|\epsilon_j$};
\draw[arr] (1.1,4.45)--(7.2,4.45) node[below] {$\delta h_1$};
\draw[arr] (4.1,1.8)--(4.1,7.15) node[left] {$\delta h_2$};
\draw[dashed,wire] (2.1,2.45) rectangle (6.1,6.45);
\draw[line width=.9pt] (2.1,4.45)--(4.1,6.45)--(6.1,4.45)--(4.1,2.45)--cycle;
\draw[arr] (4.1,4.45)--(5.1,5.45) node[above right,tinytext] {$d_1$};
\draw[arr] (4.1,4.45)--(5.1,3.45) node[below right,tinytext] {$d_2$};
\foreach \xx/\lab in {2.1/-2,6.1/2}{\draw[wire] (\xx,4.39)--(\xx,4.51);\node[tinytext,below] at (\xx,4.37){$\lab$};}
\foreach \yy/\lab in {2.45/-2,6.45/2}{\draw[wire] (4.04,\yy)--(4.16,\yy);\node[tinytext,left] at (4.0,\yy){$\lab$};}
\node[tinytext,align=center] at (3.9,1.5){$D=[d_1\ d_2]=\begin{bmatrix}1&1\\1&-1\end{bmatrix}$, $\epsilon_1=\epsilon_2=1$ (normalized).};
\node[tinytext] at (3.9,.92){The dashed box contains unattainable corners.};
\node[eq] at (12.35,8.05){$e=T^{\mathrm{record}}-T^{\mathrm{true}}$\\[5pt]$\operatorname{bias}(\widehat{\Delta h})=-l\trans\Gamma(e_t-e_b)$};
\node[eq] at (12.35,6.82){$e_t=e_b\Rightarrow\text{offset cancels}$\\[5pt]$|e_t|,|e_b|\leq\epsilon\Rightarrow|e_t-e_b|\leq2\epsilon$};
\node[eq] at (12.35,5.24){$\begin{array}{c|c|c}\text{Illustrative gain}&\text{Allocation}&\epsilon_{\mathrm{relative}}\\\hline
4\ \mathrm{mV/K}&10\ \mathrm{mV}&2.5\ \mathrm K\\0.16\ \mathrm{m}\Omega/\mathrm K&0.05\ \mathrm{m}\Omega&0.3125\ \mathrm K\end{array}$};
\node[tinytext,align=center] at (12.35,3.65){The separate SIL electrical law has 2.5 mV/K\\for the voltage contrast and zero resistance offset gain\\only because its current slopes are temperature independent.};
\node[tinytext,align=center] at (12.35,2.12){Coefficient error and causal temperature lag need their own\\discrepancy terms; junction truth is not provided by\\a common case-temperature indication.};
\draw[wire] (.1,.52)--(16.4,.52);
\node[tinytext] at (8.25,.12){$l\trans\Gamma=0$ is necessary for arbitrary offsets; the complete topology and nuisance support test remains required.};\end{tikzpicture}'''
add('05_thermal_gauge','Thermal confounding as a kernel, and calibration uncertainty as a set',b,
r'''\textbf{Figure 5.} Unbounded offsets create an exact equivalence; bounded relative offsets create a deterministic target-error set. In (b), $\alpha=2$ mV/K and the 10 mV observation are analytic settings; dashed lines mark a hypothetical $\pm1$ mV observation band, not measured noise. The diamond is a normalized linear image of an offset box. Numerical allocations belong to the illustrative law and are separated from the SIL surrogate. They are conditional budgets, not demonstrated sensor accuracies.''')

b=picture()+panel(.1,15.45,7.6,'(a) Shared-reference covariance')+panel(8.35,15.45,8.05,'(b) Correlated difference errors')
b+=r'''\node[eq] at (3.9,14.35){$\xi=(\xi_R,\xi_N,\xi_A)\trans,\quad\Cov(\xi)=\diag(v_R,v_N,v_A)$};
\node[eq] at (3.9,13.1){$\binom{\epsilon_0}{\epsilon_1}=B\xi,\qquad B=\begin{bmatrix}-1&1&0\\-1&0&1\end{bmatrix}$};
\node[eq] at (3.9,11.68){$V_\Delta=B\diag(v_R,v_N,v_A)B\trans$\\[5pt]$=\begin{bmatrix}v_N+v_R&v_R\\v_R&v_A+v_R\end{bmatrix}$};
\node[tinytext,align=center] at (3.9,10.24){Fresh reference per trial; independent record errors.\\The two contrasts within a triplet are not independent.};
\begin{axis}[mathaxis,axis equal image,at={(10.45cm,10.95cm)},anchor=south west,width=7.2cm,height=5.4cm,
 xmin=-2,xmax=2,ymin=-2,ymax=2,xtick={-1,0,1},ytick={-1,0,1},
 xlabel={Null error $\epsilon_0/\sqrt v$},ylabel={Injected error $\epsilon_1/\sqrt v$},clip=false]
\addplot[domain=0:360,samples=121,line width=.9pt]
 ({(sqrt(3)*cos(x)+sin(x))/sqrt(2)},{(sqrt(3)*cos(x)-sin(x))/sqrt(2)});
\addplot[domain=-1.6:1.6,samples=2,dashed]{x};
\addplot[domain=-1:1,samples=2,densely dotted]{-x};\end{axis}
\node[tinytext,align=center] at (12.35,9.6){$v_R=v_N=v_A=v:\ V_\Delta=v\begin{bmatrix}2&1\\1&2\end{bmatrix}$; correlation $=1/2$.};
'''
b+=panel(.1,8.95,16.3,'(c) Frozen null calibration: a marginal exchangeability statement')
b+=r'''\node[eq] at (8.25,8.1){$S_i=|\widehat{\Delta h}_{i,0}|/s_{\Delta,i},\quad n_{\mathrm{cal}}=200,\quad
k=\lceil201(1-0.05)\rceil=191,\quad q=S_{(191)}.$};
\node[eq] at (8.25,7.36){$\Pr_{\mathrm{cal,new}}\{S_{\mathrm{new}}>S_{(191)}\}\leq\frac{201-191}{201}=\frac{10}{201}<0.05$};
'''
b+=panel(.1,6.7,7.6,'(d) Gaussian reference power')+panel(8.35,6.7,8.05,'(e) Held-out software decision rates')
z=norm.ppf(.975);xs=np.linspace(0,6,61)
known=norm.cdf(-z-xs)+norm.sf(z-xs)
indep=norm.cdf(-z-xs/math.sqrt(2))+norm.sf(z-xs/math.sqrt(2))
b+=r'''\begin{axis}[mathaxis,at={(.85cm,2.6cm)},anchor=south west,width=7.15cm,height=4.0cm,
 xmin=0,xmax=6,ymin=0,ymax=1,xtick={0,2,4,6},ytick={0,.5,.8,1},
 xlabel={$\lambda=\delta/s_{\mathrm{record}}$},ylabel={Gaussian two-sided power},
 legend style={at={(.5,1.03)},anchor=south,legend columns=2,/tikz/column sep=3pt}]
'''
for curve,sty,lab in [(known,'solid','Known baseline'),(indep,'dashed','Two estimated records')]:
    coords=' '.join(f'({x:.6f},{y:.8f})' for x,y in zip(xs,curve))
    b+=r'\addplot[line width=.8pt,'+sty+'] coordinates {'+coords+r'};\addlegendentry{'+lab+'}\n'
b+=r'''\addplot[densely dotted,domain=0:6,samples=2,forget plot]{.8};
\addplot[only marks,mark=o,mark size=2pt,forget plot] coordinates {(2.801582,.8)};
\addplot[only marks,mark=*,mark size=2pt,forget plot] coordinates {(3.962034,.8)};
\end{axis}
\node[eq] at (3.9,1.08){$\pi(\tau)=\Phi(-z_{.975}-\tau)+1-\Phi(z_{.975}-\tau)$\\[3pt]
$\tau=\lambda\text{ or }\lambda/\sqrt2$ for equal independent record SE.};
'''
rows=[]
for bits in (16,12):
    for t,name in enumerate(('Voltage','Resistance')):
        d=DATA['scenarios'][str(bits)];v=d['targets'][t]
        rows.append(f"{bits}&{name}&{d['calibration_quantiles'][t]:.3f}&{100*v['calibrated_false_alarm']['rate']:.2f}&{100*v['calibrated_detection']['rate']:.2f}"+chr(92)*2)
b+=r'\node[eq] at (12.35,5.28){$\begin{array}{c|l|c|r|r}\text{ADC}&\text{Target}&q&\text{FAR [\%]}&\text{Detect [\%]}\\\hline'+'\n'+'\n'.join(rows)+r'\end{array}$};'+'\n'
b+=r'''\node[tinytext,align=center] at (12.35,3.7){200 calibration pairs and 400 evaluation triplets per ADC.\\Injected changes: 50 mV and 0.2 m$\Omega$.\\Table values are retained empirical results, not Gaussian power.};
\node[tinytext,align=center] at (12.35,2.42){16-bit median difference SE: 2.474 mV / 0.180 m$\Omega$.\\The reference contribution is included. Exact scale\\and oracle junction compensation remain model assumptions.};
\draw[wire] (.1,.55)--(16.4,.55);
\node[tinytext] at (8.25,.18){No conditional 5\% guarantee for every calibration split; no simultaneous-target or repeated-online-alarm guarantee.};\end{tikzpicture}'''
add('06_reference_decision','Baseline covariance, calibration rank and power of the actual task',b,
r'''\textbf{Figure 6.} Independent record errors induce within-triplet correlation; exchangeable null scores justify finite-rank calibration per target and ADC. The ellipse is a unit Mahalanobis contour. In (d), 80\% Gaussian power requires $\lambda\simeq2.80$ with a known baseline and $\lambda\simeq3.96$ with two equally precise independent estimates. Panel (e) reports empirical software results from quantized feedback, separately from Gaussian power. Neither calculation establishes physical calibration or online alarm validity.''')

from mathematical_visual_extensions import append_visual_extensions
append_visual_extensions(add, picture, panel)

A=np.zeros((6,8),int);A[:4,:4]=np.eye(4,dtype=int);A[4,4:6]=1;A[5,6:8]=1
K=np.zeros((8,2),int);K[4,0]=1;K[5,0]=-1;K[6,1]=1;K[7,1]=-1
h=np.array([-1,-1,1,1,-1,0,1,0]);g=np.array([1,1,1,1,1,0,1,0]);zeta=np.array([1,1,-1,-1,1,-1])
assert not np.any(A@K);assert np.array_equal(A@h,-zeta);assert np.array_equal(A@g,np.ones(6))
assert np.linalg.matrix_rank(A)==6
atlas=PRE+'\n\\begin{document}\n'
for i,f in enumerate(figures,1):
    if i>1:atlas+='\n\\newpage\n'
    atlas+=r'\PlateTitle{'+str(i)+'}{'+f['title']+'}\n'+f['body']+'\n'+r'\PlateCaption{'+f['caption']+'}\n'
    (DEST/(f['slug']+'.tikz.tex')).write_text(f['body']+'\n')
    (DEST/(f['slug']+'.tex')).write_text(PRE+'\n\\begin{document}\n'+f['body']+'\n\\end{document}\n')
atlas+='\n\\end{document}\n';(DEST/'measurement_illustrations.tex').write_text(atlas)
styles=PRE[PRE.index(r'\usetikzlibrary'):PRE.index(r'\newcommand{\PlateTitle}')].replace(r'\newcommand',r'\providecommand')
styles=styles.replace(r'\pagestyle{empty}\setlength{\parindent}{0pt}','')
styles=re.sub(r'\\DeclareMathOperator\{(\\\w+)\}\{([^}]+)\}',lambda m:r'\providecommand{'+m[1]+r'}{\operatorname{'+m[2]+'}}',styles)
(DEST/'illustration_styles.tex').write_text(styles)
(DEST/'captions.json').write_text(json.dumps({f['slug']:f['caption'] for f in figures},indent=2)+'\n')
(DEST/'provenance.json').write_text(json.dumps(dict(edition='3: mathematical geometry and visual companion plates',figures=len(figures),slugs=[f['slug'] for f in figures],
 sources=['main.tex: incidence, joint kernel, profiled SVD and shared information','thermal_identifiability.tex','measurement_baseline.tex','measurement_independent_sil.tex'],
 numerical_source='data/independent_sil_results.json',numerical_source_sha256=hashlib.sha256((ROOT/'data/independent_sil_results.json').read_bytes()).hexdigest(),
 illustrative_settings=dict(singular_values=[1,.2,0],balanced_occupancy=[.5,.5],d0=.25,thermal_scalar_gain_mV_per_K=2,thermal_scalar_observation_mV=10),
 analytical_companion_settings=dict(narrow_currents=[.8,.8,1.2,1.2],wide_currents=[0,0,2,2],
    target_tangent=[1,1,0],temperature_step_K=8,lag_gain_mV_per_K=2.5,thermal_budget_mV=5,
    lag_model='first-order analytic step; pointwise correction',per_test_alpha=.05,
    repeated_decision_curve='independent-test scenario',trajectory_protocol='proposed, not evaluated'),
 monochrome=True,hardware=False,AI_visualization_code_assistance='OpenAI Codex'),indent=2)+'\n')
(DEST/'mathematical_fixtures.json').write_text(json.dumps(dict(A=A.tolist(),K=K.tolist(),h=h.tolist(),g=g.tolist(),zeta=zeta.tolist(),
 power_lambda=xs.tolist(),power_known=known.tolist(),power_independent_equal_records=indep.tolist()),indent=2)+'\n')
print(f'{len(figures)} monochrome mathematical plates generated:',DEST)
