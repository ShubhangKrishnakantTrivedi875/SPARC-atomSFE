# Imaginary-frequency quadrature for the RPA correlation energy

Reference for the grids registered in
[`src/xc/rpa_omega_grids.py`](../../src/xc/rpa_omega_grids.py) and swept by this directory.

Each mapping gets one section with the same structure: **definition**, **constants and
reference**, **pole locations** (closed form, verified numerically), **behaviour**, and
**measured** cost. Summary and reproduction instructions are at the end.

> Rendered view: `Ctrl+Shift+V` in VS Code — the math and the figure need the preview, not the
> source.

---

## Contents

| § | mapping | status | Au nodes @1e-5 |
|---|---|---|---|
| [1](#1-rational--algebraic-mapping-of-degree-p) | rational / algebraic, degree $p$ (includes inverse-linear) | active | 300–400 |
| [2](#2-tangent-mapping) | tangent | active | 200–400 |
| [3](#3-logarithmic-mapping--withdrawn) | logarithmic | withdrawn | wrong limit |
| [4](#4-log-uniform-mapping) | log-uniform | active | 50–80 |
| [5](#5-double-exponential-mapping-de_sinh) | double-exponential | active | 30–50 |
| [6](#6-sinh-mapping-johnstonelliott) | sinh (Johnston–Elliott) | active | **20** |
| [7](#7-elliptic-sc-mapping--never-exercised) | elliptic sc | never exercised | — |
| [8](#8-gausslaguerre--withdrawn) | Gauss–Laguerre | withdrawn | overflow |
| [9](#9-aaa--no-mapping-at-all) | AAA (no mapping) | active | 30 |

---

## Notation

| symbol | meaning |
|---|---|
| $\omega\in[0,\infty)$ | imaginary frequency — the physical integration variable |
| $\tilde{\omega}\in[-1,1]$ | reference variable in which the base rule is defined |
| $w_j$ | quadrature weights in $\omega$ |
| $\tilde{w}_j$ | weights of the base rule on $[-1,1]$ |
| $c$ | scale constant of a mapping |
| $p$ | degree of the algebraic mapping |
| $n$ | number of frequency points |
| $\varepsilon_i,\ \varepsilon_a$ | occupied and partner Kohn–Sham eigenvalues |
| $\Delta_{ia}=\varepsilon_a-\varepsilon_i$ | transition energy; $\chi_0$ has poles at $\omega=\pm i\Delta_{ia}$ |
| $\Delta_{\min},\ \Delta_{\max}$ | extremes of the transition-energy band |
| $R=\Delta_{\max}/\Delta_{\min}$ | band ratio |
| $L''$ | angular coupling channel |

---

## The quantity being computed

The RPA correlation energy, summed over the angular coupling channels of the radial problem:

$$E_c^{\rm RPA}=\frac{1}{2\pi}\sum_{L''}(2L''+1)\int_0^{\infty}d\omega\;
\Big\{\ln\det\big[\mathbb{1}-\nu^{L''}\chi_0^{L''}(i\omega)\big]
+\mathrm{Tr}\big[\nu^{L''}\chi_0^{L''}(i\omega)\big]\Big\}$$

For each channel the independent-particle response is a matrix on the radial quadrature grid,
built from the Kohn–Sham orbitals $\phi$:

$$\chi_0^{L''}(r,r';i\omega)=\frac{1}{2L''+1}\sum_{i\,\in\,\rm occ}\sum_{a}
S_{ia}\begin{pmatrix}l_i & l_a & L''\\[1mm] 0 & 0 & 0\end{pmatrix}^{2}
A_{ia}\,\frac{\Delta_{ia}}{\Delta_{ia}^{2}+\omega^{2}}\;
\phi_i(r)\phi_a(r)\,\phi_i(r')\phi_a(r')$$

$$A_{ia}=f_i\,(2l_a+1)-f_a\,(2l_i+1),\qquad
S_{ia}=\begin{cases}2,& a\ \text{virtual}\\ 1,& a\ \text{occupied}\end{cases}$$

The squared Wigner 3j symbol is the angular part of the Coulomb coupling; $f$ are occupations, so
$A_{ia}$ is correct for fractional occupation and vanishes for occupied–occupied pairs with equal
$f/(2l+1)$; $S_{ia}=2$ on the occupied–virtual block accounts for the virtual–occupied block,
which is not summed separately. (Source: `_build_rpa_response_kernel` and `_build_the_constants`
in [`rpa.py`](../../src/xc/rpa.py).)

**Only one feature of this matters for the quadrature.** Every term carries the same denominator
$\Delta_{ia}^{2}+\omega^{2}$, and the orbitals, 3j symbols, $A_{ia}$, $S_{ia}$, $(2L''+1)$ and the
Coulomb kernel $\nu^{L''}$ are all **independent of $\omega$**. They set the residues, not the pole
positions. So as a function of complex $\omega$ the integrand is analytic on the real axis with
poles at

$$\omega=\pm i\Delta_{ia},\qquad \Delta_{ia}=\varepsilon_a-\varepsilon_i>0$$

and everything below depends on the spectrum only through the **set** $\{\Delta\}$ and its
extremes. For the five test elements $\Delta_{\min}$ is the HOMO–LUMO gap (0.10–0.60 Ha) and
$\Delta_{\max}$ is the top of the finite-element basis ($3.7\times10^{8}$ Ha), so $R$ reaches
$3\times10^{9}$; Au has 58 310 distinct $\Delta$.

---

## Common structure

Every grid is a **base rule** on $[-1,1]$ composed with a **mapping** $\omega(\tilde{\omega})$:

$$\int_0^\infty f(\omega)\,d\omega=\int_{-1}^{1} f\big(\omega(\tilde\omega)\big)
\frac{d\omega}{d\tilde\omega}\,d\tilde\omega
\qquad\Longrightarrow\qquad
\omega_j=\omega(\tilde\omega_j),\quad
w_j=\tilde w_j\left.\frac{d\omega}{d\tilde\omega}\right|_{\tilde\omega_j}$$

Only the mapping and its Jacobian are listed below; the weights always follow from this. Keeping
the two axes independent is deliberate — several grids described in the literature as distinct
methods are the *same* mapping at a different $c$, or the same mapping under a different base rule.

### Why the pole locus is the whole story

The base rules are exact for functions analytic in some region of the $\tilde{\omega}$-plane, and
the convergence rate is set by how far the nearest singularity of the transplanted integrand lies
from $[-1,1]$ — a Bernstein ellipse for Gauss and Clenshaw–Curtis, a strip for the trapezoidal
rule. The mapping's only job is to push the images of $\omega=\pm i\Delta$ away from $[-1,1]$,
**uniformly in $\Delta$**. A map that succeeds for one $\Delta$ but not another needs a node count
that grows with the width of the band — which is exactly the observed growth with $Z$.

<img src="mapping_pole_loci.png" alt="pole loci of chi_0 under each mapping" width="100%">

*Figure: images of the poles $\omega=\pm i\Delta$ in the $\tilde{\omega}$-plane, for $\Delta$
spanning Au's full band, coloured by $\log_{10}\Delta$. The thick black line is the integration
interval. The last panel gives distance from that interval, which is what ranks the maps. From
[`make_mapping_figures.py`](make_mapping_figures.py); every locus is the closed-form formula in
the corresponding section, verified by substituting back into the mapping (residuals
$\le2.5\times10^{-14}$).*

---

## 1. Rational / algebraic mapping of degree p

### Definition

$$\omega=c\left(\frac{1+\tilde\omega}{1-\tilde\omega}\right)^{p},
\qquad
\frac{d\omega}{d\tilde\omega}=2cp\,\frac{(1+\tilde\omega)^{p-1}}{(1-\tilde\omega)^{p+1}}$$

### Constants and reference

$c$ sets the frequency scale, $p$ concentrates nodes toward large $\omega$. The same $p=1$ map
appears in the literature at three scales:

| $c$ | source | preset |
|---|---|---|
| 2.5 | PRL **134**, 016402, supplemental material | `gl_rational_p1_a2.5` |
| 1 | Boqin MATLAB — the active path in `OneAtomFEM/basic_ingredients.py`, and the SPARC-atomSFE default | `gl_invlin` |
| 0.5 | Ren *et al.*, New J. Phys. **14**, 053020 (2012) — the FHI-aims "modified GL" grid | `gl_rational_p1_a0.5` |

**Inverse-linear is this mapping, not a separate family.** The production form
$\omega=(1-\tilde\omega)/(1+\tilde\omega)$ is $p=1$, $c=1$ under $\tilde\omega\to-\tilde\omega$,
giving the same locus mirrored. Verified: `gl_invlin` and `gl_rational_p1_a1` agree to
$8.9\times10^{-16}$ in the nodes and $3.4\times10^{-12}$ in the weights, and their measured node
counts are identical at every element.

### Pole locations

Setting $\omega=i\Delta$ gives $u=(1+\tilde\omega)/(1-\tilde\omega)=(i\Delta/c)^{1/p}$ with $p$
branches, so

$$\boxed{\;\tilde\omega_{\rm pole}
=\frac{(\Delta/c)^{1/p}\,e^{i(\pi/2+2\pi k)/p}-1}
{(\Delta/c)^{1/p}\,e^{i(\pi/2+2\pi k)/p}+1}\;},\qquad k=0,\dots,p-1$$

For $p=1$ this closes on the unit circle:

$$\tilde\omega_{\rm pole}=\frac{i\Delta-c}{i\Delta+c}
=\exp\Big\{i\big[\pi-2\arctan(\Delta/c)\big]\Big\},
\qquad |\tilde\omega_{\rm pole}|=1\ \text{exactly}$$

### Behaviour

The pole angle sweeps from $\tilde\omega=-1$ at $\Delta\ll c$, through $\tilde\omega=i$ **exactly
at $\Delta=c$** — the farthest attainable point, Bernstein radius $\rho=1+\sqrt2$ — to
$\tilde\omega=+1$ at $\Delta\gg c$. The locus is anchored to the interval endpoints, so across a
nine-decade band some pole always crowds one of them.

For $p\ge2$ the locus is a circular arc through $\pm1$ meeting the axis at angle $\pi/2p$, and its
maximum height collapses:

| $p$ | 1 | 2 | 3 | 4 |
|---|---|---|---|---|
| $\max\operatorname{Im}\tilde\omega$ | 1.00 | 0.41 | 0.27 | 0.20 |

### Measured

Nodes to reach $10^{-5}$ ($10^{-4}$ in brackets):

| preset | He | Ar | Zn | Mo | Au |
|---|---|---|---|---|---|
| `cc_rational_p1_a2.5` | 15 (10) | 80 (50) | 150 (80) | 150 (80) | 200 (150) |
| `gl_rational_p1_a2.5` | 30 (15) | 150 (80) | 200 (80) | 200 (150) | 300 (200) |
| `cc_rational_p1_a1` = `cc_invlin` | 30 (15) | 150 (80) | 200 (100) | 200 (150) | 300 (200) |
| `cc_rational_p1_a0.5` | 30 (20) | 200 (150) | 300 (150) | 300 (200) | 400 (300) |
| `gl_rational_p1_a1` = `gl_invlin` | 50 (20) | 200 (100) | 300 (150) | 300 (200) | 400 (300) |
| `gl_rational_p1_a0.5` | 50 (30) | 300 (150) | 400 (200) | 400 (300) | 400 (300) |
| $p=2,3,4$ (both rules) | \* | \* | \* | \* | \* |

Cost grows roughly linearly with $Z$, and larger $c$ helps monotonically ($0.5\to1\to2.5$ improves
at every element), consistent with the optimum sitting near the geometric centre of the band —
untested, and a cheap experiment. The $p\ge2$ grids (\*) converge smoothly to values
$5\times10^{-4}$ to $3\times10^{-3}$ from consensus and are excluded from all rankings.

---

## 2. Tangent mapping

### Definition

$$\omega=c\tan\frac{\pi(1+\tilde\omega)}{4},
\qquad
\frac{d\omega}{d\tilde\omega}=\frac{c\pi}{4}\sec^{2}\frac{\pi(1+\tilde\omega)}{4}$$

### Pole locations

With $\theta=\pi(1+\tilde\omega)/4$, $\omega=i\Delta$ requires $\tan\theta=i\Delta/c$, i.e.
$\theta=i\,\mathrm{artanh}(\Delta/c)$:

$$\boxed{\;\tilde\omega_{\rm pole}=
\begin{cases}
-1+\dfrac{4i}{\pi}\,\mathrm{artanh}\left(\Delta/c\right), & \Delta<c\\[3mm]
+1+\dfrac{4i}{\pi}\,\mathrm{artanh}\left(c/\Delta\right), & \Delta>c
\end{cases}\;}$$

### Behaviour

Poles lie on the **vertical lines $\operatorname{Re}\tilde\omega=\mp1$** at a height that vanishes
as $\Delta/c\to0$ or $\infty$ — structurally the same endpoint-anchored defect as §1, with no
compensating advantage.

### Measured

| preset | He | Ar | Zn | Mo | Au |
|---|---|---|---|---|---|
| `cc_tan_a2.5` | 20 (15) | 100 (80) | 150 (80) | 200 (100) | 200 (150) |
| `cc_tan_a1` | 30 (20) | 200 (100) | 300 (150) | 300 (200) | 400 (300) |
| `gl_tan_a2.5` | 30 (15) | 200 (80) | 300 (100) | 300 (150) | 400 (200) |
| `gl_tan_a1` | 50 (30) | 300 (150) | 300 (150) | 400 (200) | 400 (300) |

---

## 3. Logarithmic mapping — withdrawn

### Definition

$$\omega=-c\ln\frac{1-\tilde\omega}{2},
\qquad
\frac{d\omega}{d\tilde\omega}=\frac{c}{1-\tilde\omega}$$

### Pole locations

$$\boxed{\;\tilde\omega_{\rm pole}=1-2e^{-i\left(\Delta/c+2\pi k\right)},\qquad k\in\mathbb{Z}\;}$$

### Behaviour

The **circle $|\tilde\omega-1|=2$** — and the $2\pi k$ means the image *wraps around it
repeatedly* as $\Delta$ grows, so pole images accumulate densely on a curve passing through
$\tilde\omega=-1$.

### Why withdrawn

The reach is only $\omega\sim c\ln[2/(1-\tilde\omega)]$, so at $n=40$ the largest node is
$\omega\approx7$ ($c=1$) or $17$ ($c=2.5$) against a band ending at $3.7\times10^{8}$. All four
presets converged to visibly wrong limits for He — $-0.080177$ and $-0.083073$ (CC), $-0.080572$
and $-0.083203$ (GL) against the consensus $-0.0843102$, i.e. 1.3–4.9 % off — and were still
drifting at $n=400$. The mapping remains in `build_omega_grid` for exponentially decaying
integrands.

---

## 4. Log-uniform mapping

### Definition

$$\omega=e^{t},\quad
t=\ln\omega_{\min}+\tfrac12(1+\tilde\omega)L,\quad
L=\ln\frac{\omega_{\max}}{\omega_{\min}},
\qquad
\frac{d\omega}{d\tilde\omega}=\frac{\omega L}{2}$$

Nodes equispaced in $\ln\omega$, so every decade of the band is resolved equally.

### Pole locations

$\omega=i\Delta$ requires $t=\ln\Delta+i(\pi/2+2\pi k)$:

$$\boxed{\;\tilde\omega_{\rm pole}=-1+\frac{2}{L}
\left[\ln\frac{\Delta}{\omega_{\min}}+i\left(\frac{\pi}{2}+2\pi k\right)\right]\;}$$

### Behaviour

The **first mapping with a $\Delta$-independent locus**: a horizontal line at
$\operatorname{Im}\tilde\omega=\pi/L$, the same height for every pole. That is why its cost stops
growing with $Z$. But with $\omega_{\min}=10^{-9}$ and $\omega_{\max}=10^{7}$, $L=36.8$ and the
strip is only $0.085$ wide — flat, but thin. Both tails are truncated, so the range must bracket
the band with margin.

### Measured

| preset | He | Ar | Zn | Mo | Au |
|---|---|---|---|---|---|
| `gl_loguniform` | 80 (50) | 80 (50) | 80 (50) | 80 (50) | 50 (30) |
| `cc_loguniform` | 80 (50) | 80 (50) | 80 (50) | 80 (50) | 80 (50) |

The only grid that gets *cheaper* with $Z$ — a fixed decade span suits a wide band better than a
narrow one.

---

## 5. Double-exponential mapping (de_sinh)

### Definition

$$\omega=\omega_{\rm geo}\,e^{c\sinh(U\tilde\omega)},
\qquad \omega_{\rm geo}=\sqrt{\omega_{\min}\omega_{\max}},
\qquad
\frac{d\omega}{d\tilde\omega}=\omega\,cU\cosh(U\tilde\omega)$$

Two constants: $c$ sets how many decades are spanned, $U$ where the $\sinh$ saturates; both are
fixed from the requested band. The exponential covers the decades while the $\sinh$ makes **both**
tails decay doubly exponentially, so $\omega\to0$ and $\omega\to\infty$ each cost $O(1)$ nodes.

### Reference

Takahasi & Mori, Publ. RIMS Kyoto Univ. **9** (1974) 721–741. See also Trefethen & Weideman,
SIAM Review **56** (2014) 385–458, Table 14.1, where the canonical form is
$\xi=\tanh(\tfrac{\pi}{2}\sinh x)$ with rate $\exp(-\pi^2n/\log 4\pi n)$.

### Pole locations

$$\boxed{\;\tilde\omega_{\rm pole}=\frac{1}{U}\,
\mathrm{arcsinh}\left[\frac{\ln(\Delta/\omega_{\rm geo})+i(\pi/2+2\pi k)}{c}\right]\;}$$

### Behaviour

For large $|\ln\Delta|$ the $\mathrm{arcsinh}$ flattens and
$\operatorname{Im}\tilde\omega\approx\pi/(2cU|\ln\Delta|)$, so the locus **approaches the real axis
logarithmically** — measured $0.141\to0.058$ across seven decades of $\Delta$. Better than
log-uniform near the band centre, worse at its edges.

### Measured

| preset | He | Ar | Zn | Mo | Au |
|---|---|---|---|---|---|
| `mp_desinh` | 30 (20) | 30 (30) | 50 (30) | 30 (30) | 30 (30) |
| `cc_desinh` | 50 (30) | 50 (50) | 50 (50) | 50 (30) | 50 (50) |
| `gl_desinh` | 50 (30) | 50 (50) | 50 (50) | 50 (30) | 50 (30) |

The DE family is *defined* with the trapezoidal rule; pairing it with a Gaussian rule costs a
factor 1.7.

---

## 6. Sinh mapping (Johnston–Elliott)

### Definition

Natural variable $y\in[0,y_{\max}]$ with $y=\tfrac12(1+\tilde\omega)y_{\max}$ and
$y_{\max}=\mathrm{arcsinh}(\omega_{\max}/c)$:

$$\omega=c\sinh y,
\qquad
\frac{d\omega}{dy}=c\cosh y,
\qquad
\frac{d\omega}{d\tilde\omega}=\frac{y_{\max}}{2}\,c\cosh y,
\qquad c=\Delta_{\min}$$

### Constants and reference

Johnston & Elliott, *A sinh transformation for evaluating nearly singular boundary element
integrals*, Int. J. Numer. Meth. Engng **62** (2005) 564–578. For a singularity at $A+Bi$ on
$[a,b]$ their transformation is

$$x(t)=A+B\sinh\left[\frac{1-t}{2}\,\mathrm{arcsinh}\frac{a-A}{B}
+\frac{1+t}{2}\,\mathrm{arcsinh}\frac{b-A}{B}\right]$$

derived from the boundary-value problem $x'(t)=k\sqrt{(x-A)^2+B^2}$ — *keep the derivative small
where $x$ is near the singularity*. Our case is $A=0$, $B=\Delta_{\min}$,
$[a,b]=[0,\omega_{\max}]$, which reduces to the definition above **including the arcsinh range**.
So $c=\Delta_{\min}$ and $y_{\max}$ are the recipe's prescription, not tuned parameters.

### Pole locations

Using $\sinh(u\pm i\pi/2)=\pm i\cosh u$:

$$\boxed{\;y_{\rm pole}=\pm\,\mathrm{arccosh}\frac{\Delta}{c}+\frac{i\pi}{2}
\qquad\Longleftrightarrow\qquad
\tilde\omega_{\rm pole}=-1+\frac{2}{y_{\max}}
\left[\mathrm{arccosh}\frac{\Delta}{c}+\frac{i\pi}{2}\right]\;}$$

### Behaviour

$\operatorname{Im}y=\pi/2$ **exactly, for every $\Delta$** — the strip is both $\Delta$-independent
*and* the widest of any mapping here, and it is $c=\Delta_{\min}$ that places the nearest pole
precisely on its edge. In normalised units $\operatorname{Im}\tilde\omega=\pi/y_{\max}=0.189$,
versus $0.085$ for log-uniform.

Because the transplanted integrand is analytic in a strip of half-width $\pi/2$, the trapezoidal
rule converges as $\exp(-2\pi a/h)=\exp(-\pi^2 n/y_{\max})$ with $y_{\max}$ only logarithmic in the
band — hence a node count essentially independent of $Z$. This is a **strip** result, not a
Bernstein-ellipse result, so Gauss–Legendre cannot attain it.

### The reach must be limited

$\omega_{\max}$ is *not* $\Delta_{\max}$. Taking $\omega_{\max}=\Delta_{\max}=3.7\times10^{8}$
makes the rule **stall at $\sim10^{-3}$**: the solver's $\chi_0\sim1/\omega^2$ has decayed into
numerical noise there while the Jacobian weight is enormous. Measured against consensus:

| reach | He | Au |
|---|---|---|
| $6.5\times10^{4}$ | — | $3.0\times10^{-4}$ |
| $8.6\times10^{5}$ | $\mathbf{1.9\times10^{-7}}$ | $6.0\times10^{-6}$ |
| $2.5\times10^{6}$ | $1.6\times10^{-6}$ | $\mathbf{4.3\times10^{-6}}$ |
| $6.0\times10^{8}$ | $1.0\times10^{-3}$ | $7.9\times10^{-6}$ |

Hence `omega_ceiling = 1e6`, an **absolute** frequency. Tying the reach to $\Delta_{\max}$ instead
makes it $Z$-dependent and gives the widest band the narrowest grid. Trefethen & Weideman §14
document this floating-point limit as the standing drawback of exponential and double-exponential
rules. **This ceiling is the one parameter here not fixed by a published prescription.**

### Measured

| preset | He | Ar | Zn | Mo | Au |
|---|---|---|---|---|---|
| `mp_hht_elliptic` (midpoint) | **15 (15)** | **20 (15)** | **20 (15)** | **20 (15)** | **20 (15)** |
| `cc_hht_elliptic` | 15 (15) | 30 (15) | 30 (30) | 30 (20) | 30 (30) |
| `gl_hht_elliptic` | 15 (10) | 30 (20) | 30 (30) | 30 (20) | 30 (30) |

Flat in $Z$, and 20× cheaper than the production default at Au.

> **Naming.** The preset is registered as `hht_elliptic`, crediting Hale, Higham & Trefethen.
> That is wrong: their map is a Möbius transform of $\mathrm{sn}$ on a closed contour with
> $k=(\sqrt{M/m}-1)/(\sqrt{M/m}+1)$. What is owed to them is the *principles* — an elliptic-type
> map for singularities confined to a band, the $\ln 4R$ scaling, and the requirement of an
> equispaced companion rule. The map itself is Johnston & Elliott's. Rename pending.

---

## 7. Elliptic sc mapping — never exercised

### Definition

$$\omega=c\,\frac{\mathrm{sn}(y\,|\,m)}{\mathrm{cn}(y\,|\,m)},
\qquad k=\frac{\Delta_{\min}}{\Delta_{\max}},\quad m=1-k^{2},\quad c=\Delta_{\min}$$

### Pole locations

$$\boxed{\;\operatorname{Im}y_{\rm pole}=K'(m)=K(k^{2}),\qquad
\Delta\in[\Delta_{\min},\Delta_{\max}]\;}$$

verified to $\le2.5\times10^{-14}$ at $k=10^{-3}$, $0.3$ and $0.5$ — and it correctly **fails
outside the design band**: at $k=0.3$ with $\Delta=20>\Delta_{\max}=c/k$ the root jumps to a
periodic image.

### Behaviour

Since $K'(m)\to\pi/2$ as $k\to0$, this reduces to §6. With $k\approx10^{-9}$–$10^{-10}$ for all
five elements, $1-k^2$ is not representable in double precision and the implementation always
takes the $\sinh$ branch — the elliptic form is never evaluated. Its only advantage would appear
for a band narrow enough that $m$ is representable, which does not occur here.

*Historical note:* an earlier version took $K$ from `ellipk(1-k**2)`, which returns $\infty$ once
$k<10^{-8}$ and produced NaN for every element on the true band. The fix is `ellipkm1(k**2)` plus
the exact $k\to0$ closed form.

---

## 8. Gauss–Laguerre — withdrawn

No mapping: the rule is native to $[0,\infty)$ with weight $e^{-x}$, which it divides back out, so
the weights carry $e^{+x}$. The RPA integrand decays only algebraically ($\sim\omega^{-4}$ once the
$+\nu\chi_0$ counterterm cancels the leading order), so convergence was far too slow ($-0.084296$
at $n=150$ against $-0.0843102$) and for $n\ge200$ the weights overflow to non-finite $E_c$.

---

## 9. AAA — no mapping at all

### Method

Horning & Trefethen, *Quadrature formulas from rational approximations*, IMA J. Numer. Anal.
(2026), doi:10.1093/imanum/draf138. For $I=\int_\gamma f(z)w(z)\,dz$ they form the Cauchy transform
of the **weight**,

$$\mathcal{C}(s)=\frac{1}{2\pi i}\int_\gamma\frac{w(z)\,dz}{s-z},
\qquad I=\int_\Gamma f(s)\,\mathcal{C}(s)\,ds$$

AAA-approximate $r_n(s)\approx2\pi i\,\mathcal{C}(s)$ on an enclosing contour $\Gamma$, and read off

$$r_n(s)=c_\infty+\sum_k\frac{c_k}{s-z_k}
\qquad\Longrightarrow\qquad
I_n=\sum_k c_k\,f(z_k)$$

so the **poles are the nodes and the residues are the weights**. AAA never sees $f$; it enters only
by being evaluated at the nodes, and by constraining where $\Gamma$ may lie.

### How `aaa_sinh` is constructed

`_aaa_sinh_rule(n, delta_min, delta_max, omega_ceiling)` in `rpa_omega_grids.py`. The rule is
built **in the sinh variable of §6**, not in $\omega$, and the sinh map contributes only the
geometry — the nodes themselves come out of the rational approximation.

1. **Interval.** $y_{\max}=\min\!\big[\ln(4\Delta_{\max}/\Delta_{\min}),\;
   \mathrm{arcsinh}(\omega_{\max}/\Delta_{\min})\big]$, so $\gamma=[0,y_{\max}]$ with
   $\omega_{\max}=$ `omega_ceiling`. For Au, $y_{\max}=16.61$.
2. **Contour.** $\Gamma$ = stadium around $\gamma$ at half-width $h=0.9\times\pi/2$, sampled with
   300 points on each long side ($\operatorname{Im}s=\pm h$) plus 60-point semicircular caps.
   $h<\pi/2$ keeps $\Gamma$ inside the strip where $f$ is analytic — §6 put the poles at
   $\operatorname{Im}y=\pm\pi/2$.
3. **Target.** With $w=1$ on $[0,y_{\max}]$ the Cauchy transform is elementary:
   $2\pi i\,\mathcal{C}(s)=\int_0^{y_{\max}}\frac{dz}{s-z}=\log\frac{s}{s-y_{\max}}$,
   which is the $w=1$ logarithm of their §3.
4. **Fit.** AAA (Nakatsukasa–Sète–Trefethen) with `mmax = n+1`, tolerance $10^{-13}$, returning
   support points, values and barycentric weights.
5. **Poles.** Generalised eigenproblem $E\,v=\lambda B\,v$ with $B=\mathbb{1}$ except $B_{00}=0$ and
   $E=\begin{pmatrix}0&\mathbf{w}^{T}\\ \mathbf{1}&\mathrm{diag}(z_j)\end{pmatrix}$; keep the
   finite eigenvalues.
6. **Residues.** Four-point numerical formula
   $c_k=\tfrac14\sum_{m}r_n(z_k+\delta_m)\,\delta_m$ with $\delta_m=10^{-6}|z_k|e^{2\pi i m/4}$ —
   the offset must be **relative**, since the poles span decades.
7. **Filter.** Keep $-0.3<\operatorname{Re}z_k<y_{\max}+0.3$ and $|\operatorname{Im}z_k|<h$,
   discarding spurious poles.
8. **Realify.** $y_k=\operatorname{Re}z_k$, $c_k\to\operatorname{Re}c_k$. Required because the
   solver evaluates $\chi_0$ at real frequencies only; costs nothing measurable
   ($1.98\times10^{-6}\to1.62\times10^{-6}$ at degree 30).
9. **Return.** $\omega_k=\Delta_{\min}\sinh y_k$ and $w_k=c_k\,\Delta_{\min}\cosh y_k$ — the
   Jacobian is folded into the weights, since `rpa.py` forms $\sum_k w_k F(\omega_k)$ directly.

Sanity check: $\sum_k w_k\to\omega_{\max}$, which is $\int_0^{\omega_{\max}}d\omega$. Measured
$1.000\times10^{6}$ at $n=30$ for `omega_ceiling` $=10^{6}$.

### Three requirements, each established by measurement

1. **Not in $\omega$.** There $f$ is analytic only for $|\operatorname{Im}\omega|<\Delta_{\min}$
   while the arc runs to $\omega_{\max}$, so $\Gamma$ would need half-width
   $\Delta_{\min}/\omega_{\max}\approx6\times10^{-8}$ against $0.22$ in the paper's Fig. 3. The two
   sides of $\Gamma$ then differ by the branch jump $2\pi i$ — the two-branch problem their §8 says
   plain AAA handles unreliably. Confirmed: the fit residual sticks at $O(1)$ for every degree from
   20 to 150.
2. **Transform first.** Under §6's map the poles move to $\operatorname{Im}y=\pm\pi/2$ against
   $y_{\max}\approx16.6$, i.e. $\varepsilon\approx0.085$ — the same difficulty as their benchmark.
3. **Keep $w=1$, put the Jacobian in $f$.** Setting $w=c\cosh y$ fails badly (100 % error at
   $\Delta=0.5$): $\cosh$ grows like $e^{y}$, so AAA spends every node at large $\omega$. **AAA
   optimises for the weight, not for where $f$ is large.**

### Behaviour

On the synthetic $\omega^{-2}$ kernel (Au band, worst case over $\Delta=0.5\dots10^{5}$, against
the truncated-arc exact value):

| nodes | AAA | midpoint on the same map |
|---|---|---|
| ~20 | $1.4\times10^{-4}$ | $1.9\times10^{-3}$ |
| ~25 | $2.4\times10^{-5}$ | $1.2\times10^{-3}$ |
| ~30 | $\mathbf{2.0\times10^{-6}}$ | $8.4\times10^{-4}$ |
| ~40–50 | $\mathbf{1.0\times10^{-8}}$ | $3.0\times10^{-4}$ |

The midpoint rule stagnates because it over-resolves small $\Delta$ while never fixing large
$\Delta$; AAA converges uniformly. Its weights come out $\approx y_{\max}/n$ in the interior with
tapered ends — it rediscovers the trapezoidal rule *plus* the end corrections.

### Measured

| preset | He | Ar | Zn | Mo | Au |
|---|---|---|---|---|---|
| `aaa_sinh` | 20 (20) | 20 (20) | 20 (15) | 30 (20) | 30 (15) |

Second only to §6 on the real calculation, despite winning decisively on the synthetic kernel —
the real integrand decays as $\omega^{-4}$, which is far more forgiving of large-$\Delta$ error
than the $\omega^{-2}$ test.

The AAA implementation is Nakatsukasa, Sète & Trefethen, SIAM J. Sci. Comput. **40** (2018) A1494,
validated against the §3 benchmark of the quadrature paper ($1.7\times10^{-5}$ at degree 20,
against the $1.6\times10^{-4}$ they report).

---

## Base rules

| rule | notes |
|---|---|
| `gauss_legendre` | interior nodes; optimal for functions analytic in a Bernstein ellipse |
| `clenshaw_curtis` | keeps one endpoint, giving an explicit $\omega=0$ node. Endpoint weights $1/(m^2-1)$ for even $m$, $1/m^2$ for odd. One endpoint must be dropped since every mapping sends it to $\omega=\infty$ — `drop='low'` for inverse-linear, `'high'` otherwise |
| `midpoint` | equispaced, $\tilde w=2/n$. Required by §6–§7: the transplanted integrand is periodic in a strip, where the trapezoidal rule converges as $\exp(-2\pi a/h)$. Midpoint rather than closed trapezoid because $y=y_{\max}$ maps to $\omega=\infty$ |
| `gauss_laguerre` | native to $[0,\infty)$, ignores the mapping — withdrawn |
| `aaa` | generates its own nodes and weights (§9) |

Measured effect of the base rule at fixed mapping: 20 vs 30 nodes (§6), 30 vs 50 (§5) — a
consistent 1.5–1.7× for the equispaced rule.

---

## Summary

Distance of the pole locus from $[-1,1]$ predicts the measured ordering exactly:

| mapping | $\operatorname{Im}\tilde\omega$ at the poles | $\Delta$-independent? | Au nodes @1e-5 |
|---|---|---|---|
| sinh (§6) | $\pi/y_{\max}=0.189$ | **yes** | **20** |
| AAA (§9) | — (no map) | — | 30 |
| de_sinh (§5) | $0.141\to0.058$ | no, shrinks | 30–50 |
| log-uniform (§4) | $\pi/L=0.085$ | yes | 50–80 |
| rational $p=1$ (§1) | 1 at $\Delta=c$, $\to0$ at both extremes | no, collapses | 300–400 |
| tan (§2) | $\to0$ at both extremes | no | 200–400 |
| rational $p\ge2$ (§1) | $\le0.41$, $\to0$ | no | wrong limit |
| log (§3) | wraps, accumulates | no | wrong limit |

### Consensus energies

Median over the 21 trusted presets, each self-referenced against its own $n=400$:

| element | $Z$ | $E_c^{\rm RPA}$ (Ha) |
|---|---|---|
| He | 2 | $-0.084310216825$ |
| Ar | 18 | $-1.109310300010$ |
| Zn | 30 | $-2.294550802732$ |
| Mo | 42 | $-3.072565501340$ |
| Au | 79 | $-6.803723277190$ |

**Caveat.** The trusted presets for Au span $\sim4\times10^{-5}$, while §6 and §9 already reach
$3\times10^{-8}$. Node counts at $10^{-5}$ and $10^{-4}$ are therefore solid, but any quoted
accuracy below $\sim10^{-5}$ for the heavier elements is limited by the *reference*, not by the
grid. Sharpening it needs a high-$n$ run on two or three independent grids.

---

## Reproducing

Fixed parameters: $N_{fe}=10$, polynomial order 20, quadrature order 45, $L_{\max}=20$,
`xc = RPA@DFT` on GGA-PBE orbitals, `enable_parallelization = True`.

```bash
cd tests/rpa_omega_quadrature

python save_spectra.py                       # stage 1: GGA-PBE spectra -> gga_pbe_spectra/
python _run_one.py Au 20 mp_hht_elliptic     # one case
python make_mapping_figures.py               # this document's figure

# full sweep: one SLURM job per (element, preset), each running the whole n ladder
for EL in He Ar Zn Mo Au; do
  ELEMENT=$EL PRESETS="aaa_sinh mp_hht_elliptic gl_hht_elliptic mp_desinh gl_desinh gl_invlin gl_rational_p1_a2.5" \
    sbatch --array=0-6 submit_sweep.sbatch
done

for EL in He Ar Zn Mo Au; do python merge_results.py $EL; done
python summarise_all.py                      # all active presets -> convergence_curves.pdf
python summarise_all.py --presets aaa_sinh mp_hht_elliptic mp_desinh gl_invlin
```

The $n$ ladder is 5, 10, 15, 20, 30, 50, 80, 100, 150, 200, 300, 400. Each preset is scored
against its own $n=400$, with a cross-preset consensus check to catch convergence to a wrong limit.

---

## References

1. N. Johnston and D. Elliott, *A sinh transformation for evaluating nearly singular boundary
   element integrals*, Int. J. Numer. Meth. Engng **62** (2005) 564–578.
2. L. N. Trefethen and J. A. C. Weideman, *The exponentially convergent trapezoidal rule*,
   SIAM Review **56** (2014) 385–458.
3. H. Takahasi and M. Mori, *Double exponential formulas for numerical integration*,
   Publ. RIMS Kyoto Univ. **9** (1974) 721–741.
4. A. Horning and L. N. Trefethen, *Quadrature formulas from rational approximations*,
   IMA J. Numer. Anal. (2026), doi:10.1093/imanum/draf138.
5. Y. Nakatsukasa, O. Sète and L. N. Trefethen, *The AAA algorithm for rational approximation*,
   SIAM J. Sci. Comput. **40** (2018) A1494–A1522.
6. N. Hale, N. J. Higham and L. N. Trefethen, *Computing $A^\alpha$, $\log(A)$, and related matrix
   functions by contour integrals*, SIAM J. Numer. Anal. **46** (2008) 2505–2523.
7. N. Hale and L. N. Trefethen, *New quadrature formulas from conformal maps*,
   SIAM J. Numer. Anal. **46** (2008) 930–948.
8. X. Ren *et al.*, *Random-phase approximation and its applications in computational chemistry
   and materials science*, New J. Phys. **14** (2012) 053020.
9. PRL **134**, 016402, supplemental material (`scrpa4_SM.pdf`).

*Reference 1 has not been read directly; the transformation and its BVP derivation are taken from
the restatement in arXiv:2210.09954. Confirm against the original before publication.*
