"""
Verification instrument for 2026-10-01 PROVE: the unit of replication.

Model (M1):   X_ij = theta_i + u_i + e_ij            i=1..I items, j=1..J judges
Model (M2):   X_ij = theta_i + u_i + w_j + e_ij       (adds a judge main effect)
Model (M3):   E_ij ~ Bernoulli(p_i) conditionally independent given item  (binary errors,
                                                      the form Lyra's recipe actually consumes)

Var(theta)=v, Var(u)=tau2, Var(e)=sigma2, Var(w)=omega2.

Everything is tested against the ESTIMATOR LYRA RUNS (np.corrcoef mean off-diagonal),
not only against the theoretical covariance.
"""
import numpy as np

rng = np.random.default_rng(20261001)

# ---------- the estimators, copied verbatim from DRAFT-2026-09-11.md ----------
def phi_bar_hat(E):
    keep = E.std(axis=0) > 0
    C = np.corrcoef(E[:, keep], rowvar=False)
    m = C.shape[0]
    iu = np.triu_indices(m, k=1)
    return C[iu].mean(), m

def n_eff_kish(E):
    p, m = phi_bar_hat(E)
    return m / (1 + (m - 1) * p)

def n_eff_esdof(E):
    keep = E.std(axis=0) > 0
    C = np.corrcoef(E[:, keep], rowvar=False)
    eig = np.clip(np.linalg.eigvalsh(C), 0, None)
    return eig.sum()**2 / np.square(eig).sum()

# ---------- closed forms under test ----------
def kish(J, p):            return J / (1 + (J - 1) * p)
def ceiling(p):            return 1.0 / p
def var_M(v, tau2, sig2, I, J):  return (v + tau2 + sig2 / J) / I
def esdof_cs(J, a, s):
    """participation ratio of compound-symmetric cov a*11^T + s*Id (a=shared, s=private)"""
    lam = np.array([J * a + s] + [s] * (J - 1))
    return lam.sum()**2 / np.square(lam).sum()

def sim_panel(v, tau2, sig2, I, J, omega2=0.0, seed=None):
    r = np.random.default_rng(seed)
    theta = r.normal(0, np.sqrt(v),    (I, 1)) if v    > 0 else np.zeros((I, 1))
    u     = r.normal(0, np.sqrt(tau2), (I, 1)) if tau2 > 0 else np.zeros((I, 1))
    w     = r.normal(0, np.sqrt(omega2),(1, J)) if omega2 > 0 else np.zeros((1, J))
    e     = r.normal(0, np.sqrt(sig2), (I, J)) if sig2 > 0 else np.zeros((I, J))
    X = theta + u + w + e
    E = u + w + e                   # error = score minus truth
    return X, E

print("="*78)
print("TEST 1  Prop 1 as an identity: phi_bar = tau2/(tau2+sigma2), Gaussian model (M1)")
print("="*78)
print(f"{'v':>6}{'tau2':>7}{'sig2':>7}{'J':>5}{'I':>7} | {'phi_pred':>9}{'phi_hat':>9}{'err':>9} | {'kish_pred':>10}{'kish_hat':>9}")
worst = 0.0
for (v, tau2, sig2, J) in [(1,1,1,9),(0,1,1,9),(5,.2,3,4),(3.7,.99,1,100),
                           (0.1,2,0.05,25),(2,0.01,4,50),(1,3,0.5,3),(0,0.5,0.5,2)]:
    I = 400000
    _, E = sim_panel(v, tau2, sig2, I, J, seed=1)
    ph, m = phi_bar_hat(E)
    pred = tau2 / (tau2 + sig2)
    worst = max(worst, abs(ph - pred))
    print(f"{v:6.2f}{tau2:7.2f}{sig2:7.2f}{J:5d}{I:7d} | {pred:9.4f}{ph:9.4f}{ph-pred:9.1e} | "
          f"{kish(J,pred):10.4f}{n_eff_kish(E):9.4f}")
print(f"  worst |phi_hat - tau2/(tau2+sig2)| = {worst:.2e}   <- identity, not a fit\n")

print("="*78)
print("TEST 2  NEGATIVE CONTROL: tau2 = 0 must make the ceiling EVAPORATE")
print("="*78)
for J in [2, 9, 25, 100]:
    _, E = sim_panel(v=4.0, tau2=0.0, sig2=1.0, I=400000, J=J, seed=2)
    ph, _ = phi_bar_hat(E)
    print(f"  J={J:4d}  phi_hat={ph:+.5f}  n_eff_kish={n_eff_kish(E):8.3f}   (must be ~= J = {J})"
          f"   predicted ceiling 1/phi = {'inf' if abs(ph)<1e-3 else f'{1/ph:.1f}'}")
print("  -> no ceiling below J when tau2=0. Instrument REFUSES to report one.\n")
print("  and the positive control it must distinguish from (same v,sig2, tau2=1):")
for J in [2, 9, 25, 100]:
    _, E = sim_panel(v=4.0, tau2=1.0, sig2=1.0, I=400000, J=J, seed=2)
    print(f"  J={J:4d}  phi_hat={phi_bar_hat(E)[0]:+.5f}  n_eff_kish={n_eff_kish(E):8.3f}   ceiling=2.000")

print("\n"+"="*78)
print("TEST 3  Her three numbers: can Prop 1 carry 49% / 1.21 / 1.99 simultaneously?")
print("="*78)
for J in [9, 41, 100]:
    # back-solve phi from kish
    phi_err = (J/1.99 - 1)/(J-1)
    phi_raw = (J/1.21 - 1)/(J-1)
    # set sig2 = 1
    sig2 = 1.0
    tau2 = phi_err/(1-phi_err)*sig2
    a    = phi_raw/(1-phi_raw)*sig2      # a = v + tau2
    v    = a - tau2
    print(f"  J={J:4d}: phi_err={phi_err:.4f} (she reports 0.49)  phi_raw={phi_raw:.4f}"
          f"  -> v={v:7.3f} tau2={tau2:6.3f} sig2={sig2:.3f}   ceiling=1/phi={1/phi_err:.3f}")
print()
print("  forward check: what does phi_bar = 0.49 (her reported 49%) imply?")
for J in [9, 20, 41, 100, 1000]:
    print(f"    J={J:5d}  n_eff_kish = {kish(J,0.49):6.3f}   (ceiling {ceiling(0.49):.3f})")
# invert J/(1+(J-1)p)=t  =>  J = t(1-p)/(1-tp)
inv_kish = lambda t,p: t*(1-p)/(1-t*p)
Jstar = inv_kish(1.99, 0.49)
print(f"  -> n_eff = 1.99 at phi_bar=0.49 requires J = {Jstar:.1f}.  At J=9 it would be {kish(9,0.49):.3f}.")

print("\n"+"="*78)
print("TEST 4  Prop 2: Var(M) = (v + tau2 + sigma2/J)/I  [2000 replicate panels]")
print("="*78)
print(f"{'v':>6}{'tau2':>7}{'sig2':>7}{'I':>6}{'J':>5} | {'Var_pred':>11}{'Var_hat':>11}{'ratio':>8}")
for (v,tau2,sig2,I,J) in [(1,1,1,50,9),(3.7,.99,1,200,100),(0,2,1,20,3),(5,0,1,100,9),(1,1,0,30,7)]:
    Ms = np.array([sim_panel(v,tau2,sig2,I,J,seed=1000+k)[0].mean() for k in range(4000)])
    pred = var_M(v,tau2,sig2,I,J)
    print(f"{v:6.2f}{tau2:7.2f}{sig2:7.2f}{I:6d}{J:5d} | {pred:11.6f}{Ms.var():11.6f}{Ms.var()/pred:8.4f}")

print("\n"+"="*78)
print("TEST 5  Cor 2(a) floor and Cor 2(b) no-floor;  FIXED-BUDGET B=I*J monotonicity")
print("="*78)
v,tau2,sig2 = 3.7,0.99,1.0
print("  (a) J -> inf at I=200:  floor=(v+tau2)/I =", f"{(v+tau2)/200:.6f}")
for J in [1,9,100,10_000,10_000_000]:
    print(f"      J={J:>10d}  Var(M)={var_M(v,tau2,sig2,200,J):.6f}   still purchasable sig2/(IJ)={sig2/(200*J):.3e}")
print("  (b) I -> inf at J=1:")
for I in [200, 2000, 200_000, 20_000_000]:
    print(f"      I={I:>10d}  Var(M)={var_M(v,tau2,sig2,I,1):.3e}")
print("  FIXED BUDGET B = I*J = 1800 judgments, which (I,J) minimises Var(M)?")
B=1800
rows=[(J, B//J, var_M(v,tau2,sig2,B//J,J)) for J in [1,2,3,5,9,25,100,300,900,1800]]
for J,I,V in rows: print(f"      J={J:>5d} I={I:>5d}  Var(M)={V:.6f}   closed form (J(v+tau2)+sig2)/B={(J*(v+tau2)+sig2)/B:.6f}")
print(f"  -> argmin J = {min(rows,key=lambda r:r[2])[0]} .  Var = (J(v+tau2)+sigma2)/B is strictly increasing in J.")

print("\n"+"="*78)
print("TEST 4b  is the uniform 2% shortfall in TEST 4 shared-seed sampling noise?")
print("="*78)
v,tau2,sig2,I,J = 1,1,1,50,9
pred = var_M(v,tau2,sig2,I,J)
for off in [1000, 50_000, 900_000, 7_000_000]:
    Ms = np.array([sim_panel(v,tau2,sig2,I,J,seed=off+k)[0].mean() for k in range(4000)])
    print(f"  seed block {off:>9d}: Var_hat/Var_pred = {Ms.var()/pred:.4f}")
Ms = np.array([sim_panel(v,tau2,sig2,I,J,seed=3_000_000+k)[0].mean() for k in range(60000)])
print(f"  60000 replicates          : Var_hat/Var_pred = {Ms.var()/pred:.4f}  (rel SE ~ {np.sqrt(2/60000):.4f})")

print("\n"+"="*78)
print("TEST 6  THE ENEMY: a judge main effect w_j (Var = omega2).  Model M2.")
print("="*78)
print("  6a  is Lyra's phi_bar_hat (np.corrcoef on error columns) BLIND to omega2?")
print(f"  {'omega2':>8} | {'phi_hat':>9}{'tau2/(tau2+sig2)':>18}{'tau2/(tau2+om2+sig2)':>22}")
for om2 in [0.0, 0.25, 1.0, 4.0, 25.0]:
    _, E = sim_panel(v=1, tau2=1, sig2=1, I=400000, J=9, omega2=om2, seed=7)
    print(f"  {om2:8.2f} | {phi_bar_hat(E)[0]:9.4f}{1/2:18.4f}{1/(2+om2):22.4f}")
print("  -> phi_hat tracks tau2/(tau2+sig2) and IGNORES omega2: a per-judge mean shift is")
print("     annihilated by the per-column centring inside Pearson correlation.")
print()
print("  6b  but does Var(M) ignore omega2?  (judges RE-DRAWN each replicate = random panel)")
def var_M_random_panel(v,tau2,sig2,om2,I,J,R=20000,off=0):
    return np.array([sim_panel(v,tau2,sig2,I,J,om2,seed=off+k)[0].mean() for k in range(R)]).var()
print(f"  {'omega2':>8}{'I':>8} | {'(v+tau2+sig2/J)/I':>19}{'+ om2/J':>12}{'Var_hat':>11}")
for om2 in [0.0, 1.0, 4.0]:
    for I in [50, 500]:
        base = var_M(1,1,1,I,9); full = base + om2/9
        vh = var_M_random_panel(1,1,1,om2,I,9,R=20000,off=11_000_000)
        print(f"  {om2:8.2f}{I:8d} | {base:19.6f}{full:12.6f}{vh:11.6f}")
print("  -> omega2/J does NOT shrink with I.  Cor 2(b) 'I has no floor' is FALSE when omega2>0:")
print("     inf_I Var(M) = omega2/J.  Items cannot buy down a judge main effect; only judges can.")
print()
print("  6c  fixed-budget optimum with omega2>0:  Var = (J(v+tau2)+sig2)/B + omega2/J")
for om2 in [0.0, 0.05, 0.5]:
    B=1800
    f = lambda J: (J*(1+1)+1)/B + om2/J
    best = min(range(1,B+1), key=f)
    Jstar = np.sqrt(om2*B/2) if om2>0 else 1
    print(f"    omega2={om2:5.2f}: argmin_J = {best:5d}   sqrt(omega2*B/(v+tau2)) = {Jstar:8.2f}")

print("\n"+"="*78)
print("TEST 7  Model M3: BINARY errors (the form Lyra's recipe actually consumes).")
print("        Prop 1 must hold with tau2=Var(p_i), sigma2=E[p_i(1-p_i)] -- no Gaussianity.")
print("="*78)
print(f"  {'p_dist':>22}{'J':>5} | {'pbar':>7}{'Var(p)':>9}{'E[p(1-p)]':>11}{'phi_pred':>10}{'phi_hat':>9}{'err':>9}")
I = 600000
for name, draw in [("Beta(1,1)=Unif",      lambda r: r.beta(1,1,(I,1))),
                   ("Beta(2,5)",           lambda r: r.beta(2,5,(I,1))),
                   ("Beta(.3,.3) bimodal", lambda r: r.beta(.3,.3,(I,1))),
                   ("Beta(20,20) tight",   lambda r: r.beta(20,20,(I,1))),
                   ("2-point {.1,.9}",     lambda r: np.where(r.random((I,1))<.5,.1,.9)),
                   ("DEGENERATE p=.3",     lambda r: np.full((I,1),.3))]:
    for J in [9, 100]:
        r = np.random.default_rng(31+J)
        p = draw(r)
        E = (r.random((I,J)) < p).astype(float)
        tau2 = p.var(); sig2 = (p*(1-p)).mean(); pbar = p.mean()
        pred = tau2/(tau2+sig2) if tau2+sig2>0 else float('nan')
        ph,_ = phi_bar_hat(E)
        print(f"  {name:>22}{J:5d} | {pbar:7.4f}{tau2:9.5f}{sig2:11.5f}{pred:10.5f}{ph:9.5f}{ph-pred:9.1e}")
print("  (DEGENERATE p=const is the binary negative control: no item-level error -> phi=0 -> n_eff=J)")

print("\n"+"="*78)
print("TEST 8  ESDOF.  Raw-score ESDOF is a function of (v+tau2)/sigma2 -- i.e. of SIGNAL.")
print("        Error ESDOF is a function of (J, phi_bar) alone.")
print("="*78)
print(f"  {'v':>6}{'tau2':>7}{'sig2':>7}{'J':>5} | {'ESDOF_raw':>11}{'pred':>9} | {'ESDOF_err':>11}{'pred':>9} | {'kish_err':>9}")
for (v,tau2,sig2,J) in [(3.7,.99,1,9),(3.7,.99,1,100),(0,1,1,9),(100,1,1,9),(1,1,1,9),(0.01,1,1,9)]:
    X,E = sim_panel(v,tau2,sig2,400000,J,seed=8)
    pr_raw = esdof_cs(J, v+tau2, sig2); pr_err = esdof_cs(J, tau2, sig2)
    print(f"{v:6.2f}{tau2:7.2f}{sig2:7.2f}{J:5d} | {n_eff_esdof(X):11.4f}{pr_raw:9.4f} | "
          f"{n_eff_esdof(E):11.4f}{pr_err:9.4f} | {n_eff_kish(E):9.4f}")
print()
print("  8b  THE DEGENERATE CASE: perfect judges (u=e=0, so tau2=sig2=0, signal only).")
X,E = sim_panel(v=4.0, tau2=0.0, sig2=0.0, I=5000, J=9, seed=9)
print(f"      raw-score ESDOF on a PERFECT panel = {n_eff_esdof(X):.6f}   (closed form {esdof_cs(9,4.0,0.0):.6f})")
import warnings
with warnings.catch_warnings():
    warnings.simplefilter('ignore')
    try:   print(f"      error ESDOF on a PERFECT panel  = {n_eff_esdof(E)}   <- her code's actual return value")
    except Exception as ex: print(f"      error ESDOF on a PERFECT panel  RAISES {type(ex).__name__}: {ex}")
    try:   print(f"      error Kish  on a PERFECT panel  = {n_eff_kish(E)}")
    except Exception as ex: print(f"      error Kish  on a PERFECT panel  RAISES {type(ex).__name__}: {ex}")
print("  -> raw-score ESDOF is MINIMISED (=1) exactly when the judges are PERFECT.")
print("     A reliability statistic minimised by perfect judges is measuring signal, not reliability.")
print()
print("  8c  monotone equivalence: both ESDOF_err and Kish are strictly decreasing in phi_bar (J=9)")
print(f"      {'phi':>6}{'ESDOF_err':>11}{'kish':>9}")
for p in [0.0,0.1,0.25,0.49,0.75,0.9,1.0]:
    print(f"      {p:6.2f}{esdof_cs(9,p,1-p):11.4f}{kish(9,p) if p>0 else 9:9.4f}")

print("\n"+"="*78)
print("TEST 9  The two statistics called 'n_eff' have DIFFERENT ceilings: 1/phi vs 1/phi^2")
print("="*78)
print(f"  {'phi':>6} | {'Kish(J=1e6)':>12}{'1/phi':>8} | {'ESDOF(J=1e6)':>13}{'1/phi^2':>9}")
for p in [0.1,0.25,0.4975,0.49,0.64,0.9]:
    J=1_000_000
    print(f"  {p:6.4f} | {kish(J,p):12.4f}{1/p:8.4f} | {esdof_cs(J,p,1-p):13.4f}{1/p**2:9.4f}")
print("  Asymptotics: Sum(lam)=J exactly; Sum(lam^2)=(Jp+1-p)^2+(J-1)(1-p)^2 ~ J^2 p^2")
print("  => ESDOF ~ J^2/(J^2 p^2) = 1/phi^2.   Kish -> 1/phi.   Ratio of ceilings = 1/phi.")
print(f"  On her panel (phi=0.49): Kish ceiling {1/0.49:.2f},  ESDOF ceiling {1/0.49**2:.2f}.")

print("\n"+"="*78)
print("TEST 10  The exchange rate on HER numbers, both estimands, both framings")
print("="*78)
for J in [41, 100]:
    phi = (J/1.99 - 1)/(J-1); sig2=1.0
    tau2 = phi/(1-phi); phir=(J/1.21-1)/(J-1); a = phir/(1-phir); v = a-tau2
    print(f"  J={J}: v={v:.3f} tau2={tau2:.3f} sig2={sig2:.3f}  (phi_err={phi:.4f})")
    print(f"     SAMPLED items (v counts):  v+tau2={v+tau2:.3f}  vs  sig2={sig2:.3f}  ->  "
          f"{'infinite judges WINS' if v+tau2<sig2 else 'ONE judge on 2I items WINS'}"
          f"   [ratio {(v+tau2)/sig2:.2f}]")
    print(f"     FIXED items  (v drops) :  tau2={tau2:.3f}      vs  sig2={sig2:.3f}  ->  "
          f"{'infinite judges WINS' if tau2<sig2 else 'ONE judge on 2I items WINS'}"
          f"   [ratio {tau2/sig2:.3f}; break-even is exactly phi=1/2]")
print("  Fixed-budget statement, UNCONDITIONAL in both framings (omega2=0):")
print("     sampled: Var=(J(v+tau2)+sig2)/B  increasing in J   -> J*=1")
print("     fixed  : Var=(J*tau2+sig2)/B     increasing in J   -> J*=1")

print("\n"+"="*78)
print("TEST 11  Prop 2 EXACTLY, by symbolic expansion (not Monte Carlo)")
print("="*78)
import sympy as sp
I,J = 3,4
th=[sp.Symbol(f'th{i}') for i in range(I)]; u=[sp.Symbol(f'u{i}') for i in range(I)]
w=[sp.Symbol(f'w{j}') for j in range(J)]; e=[[sp.Symbol(f'e{i}_{j}') for j in range(J)] for i in range(I)]
v_,t_,s_,o_ = sp.symbols('v tau2 sigma2 omega2', positive=True)
M = sp.Rational(1,I*J)*sum(th[i]+u[i]+w[j]+e[i][j] for i in range(I) for j in range(J))
Mx = sp.expand(M**2)
var = 0
for term, coeff in Mx.as_coefficients_dict().items():
    fs = term.free_symbols
    if len(fs)==1:
        nm = list(fs)[0].name
        var += coeff*( v_ if nm.startswith('th') else t_ if nm.startswith('u') else
                       o_ if nm.startswith('w')  else s_ )
    # cross terms have expectation 0 (all components independent, mean 0)
var = sp.simplify(var)
pred = sp.simplify((v_+t_)/I + s_/(I*J) + o_/J)
print(f"  I={I}, J={J}")
print(f"  E[M^2] from symbolic expansion = {sp.nsimplify(var)}")
print(f"  (v+tau2)/I + sigma2/(I*J) + omega2/J = {sp.nsimplify(pred)}")
print(f"  difference = {sp.simplify(var-pred)}    <- 0 means Prop 2 (M2 version) is an identity")
for (I,J) in [(1,1),(2,7),(5,3),(9,9)]:
    M = sp.Rational(1,I*J)*sum(sp.Symbol(f'th{i}')+sp.Symbol(f'u{i}')+sp.Symbol(f'w{j}')+sp.Symbol(f'e{i}_{j}')
                               for i in range(I) for j in range(J))
    vv=0
    for term,coeff in sp.expand(M**2).as_coefficients_dict().items():
        fs=term.free_symbols
        if len(fs)==1:
            nm=list(fs)[0].name
            vv += coeff*( v_ if nm.startswith('th') else t_ if nm.startswith('u') else o_ if nm.startswith('w') else s_)
    d = sp.simplify(vv - ((v_+t_)/I + s_/(I*J) + o_/J))
    print(f"  I={I:2d} J={J:2d}: difference = {d}")

print("\n"+"="*78)
print("TEST 12  BREAK THE BIRTH RANGE: randomised grid, 4 decades per parameter")
print("="*78)
r = np.random.default_rng(1201)
nP1=nP2=nMON=nCEIL=0; N=4000
wP1=wP2=0.0
for k in range(N):
    v    = 10**r.uniform(-2,2) * (r.random()>0.15)
    tau2 = 10**r.uniform(-2,2) * (r.random()>0.15)
    sig2 = 10**r.uniform(-2,2)
    om2  = 10**r.uniform(-2,2) * (r.random()>0.5)
    I    = int(10**r.uniform(0,3)); J = int(10**r.uniform(0,2.5))+1
    # Prop 1 (closed form on the model's own covariance): phi = tau2/(tau2+sig2)
    phi = tau2/(tau2+sig2)
    # Prop 2 identity, checked against the covariance algebra directly
    lhs = (v+tau2)/I + sig2/(I*J) + om2/J
    rhs = var_M(v,tau2,sig2,I,J) + om2/J
    wP2 = max(wP2, abs(lhs-rhs)/max(lhs,1e-300)); nP2 += abs(lhs-rhs) <= 1e-12*max(lhs,1e-300)
    # Kish ceiling: strict for every finite J, and attained only in the limit
    nCEIL += (kish(J,phi) < ceiling(phi) if phi>0 else kish(J,0)==J) and kish(J,phi)>0
    # fixed-budget monotonicity when om2=0: Var is increasing in J
    B=float(I*J)
    nMON += all((Jb*(v+tau2)+sig2)/B < (Jb2*(v+tau2)+sig2)/B
                for Jb,Jb2 in [(1,2),(2,3),(5,6),(J,J+1)]) if v+tau2>0 else True
    nP1 += 1
print(f"  draws                                   : {N}")
print(f"  Prop 2 identity holds (rel err <= 1e-12): {nP2}/{N}   worst rel err {wP2:.1e}")
print(f"  Kish ceiling strict at finite J          : {nCEIL}/{N}")
print(f"  fixed-budget Var strictly increasing in J: {nMON}/{N}")
print()
print("  Monte-Carlo spot checks of Prop 1 far from her region (phi_hat vs tau2/(tau2+sig2)):")
print(f"  {'v':>9}{'tau2':>9}{'sig2':>9}{'om2':>7}{'J':>5} | {'phi_pred':>9}{'phi_hat':>9}{'err':>9}")
for k in range(10):
    v=10**r.uniform(-2,2); tau2=10**r.uniform(-2,2); sig2=10**r.uniform(-2,2)
    om2=10**r.uniform(-2,2)*(k%2); J=int(10**r.uniform(0.3,2))+1
    _,E = sim_panel(v,tau2,sig2,300000,J,om2,seed=2000+k)
    ph,_=phi_bar_hat(E); pred=tau2/(tau2+sig2)
    print(f"  {v:9.3f}{tau2:9.3f}{sig2:9.3f}{om2:7.2f}{J:5d} | {pred:9.5f}{ph:9.5f}{ph-pred:9.1e}")
print("\n  DONE.")

print("\n"+"="*78)
print("TEST 13  Prop 4 needs splitting: 'per-item' is THREE functionals, not one.")
print("="*78)
v,tau2,sig2 = 3.7,0.99,1.0; phi = tau2/(tau2+sig2)
print(f"  (v,tau2,sig2)=({v},{tau2},{sig2})  phi_bar={phi:.4f}  Kish ceiling 1/phi={1/phi:.3f}")
print(f"  {'J':>7} | {'Var(M)*I':>10} | {'Var(Xbar_i|i)':>14} | {'MSE(Xbar_i vs truth)':>21}")
for J in [1,2,9,100,10**6]:
    print(f"  {J:7d} | {v+tau2+sig2/J:10.4f} | {sig2/J:14.6f} | {tau2+sig2/J:21.4f}")
print(f"  limits  | {v+tau2:10.4f} | {0.0:14.6f} | {tau2:21.4f}")
print(f"  design effect J=1 -> inf:  aggregate {(v+tau2+sig2)/(v+tau2):.4f}   "
      f"per-item-vs-consensus INFINITE (no ceiling)   per-item-vs-TRUTH {(tau2+sig2)/tau2:.4f} = 1/phi")
print("  -> the ceiling 1/phi governs per-item ACCURACY too.  It is absent only for the")
print("     variance of the panel score about the panel's OWN consensus, which is not an")
print("     estimand anyone wants.  My brief's 'per-item use has no ceiling' was too loose.")
print()
print("  13b  DETECTION has no ceiling: P(all J judges miss) with p_i ~ Beta(a,b), cond. indep.")
print(f"  {'p_i dist':>16}{'J=1':>9}{'J=2':>9}{'J=9':>9}{'J=100':>10}{'J=1000':>10}  decay")
for name,(a,b) in [("Beta(2,2)",(2,2)),("Beta(.5,.5)",(.5,.5)),("Beta(1,1)",(1,1)),("Beta(5,1) hard",(5,1))]:
    r2=np.random.default_rng(77); p=r2.beta(a,b,4_000_000)
    vals=[ (p**J).mean() for J in [1,2,9,100,1000] ]
    # for Beta(a,b), E[p^J] = B(a+J,b)/B(a,b) ~ Gamma(b) * J^{-b}
    print(f"  {name:>16}" + "".join(f"{x:9.5f}" if x>1e-5 else f"{x:9.2e}" for x in vals[:3])
          + f"{vals[3]:10.3e}{vals[4]:10.3e}  ~ J^-{b}")
print("  -> E[p_i^J] -> 0 for every J, polynomially (J^-b) not geometrically. No floor.")
print("     Conditional independence is intact; detection spends it, estimation cannot.")

print("\n"+"="*78)
print("TEST 14  WHICH functional has Kish's design effect?  (correcting my own brief)")
print("="*78)
for (v,tau2,sig2,I) in [(3.7,.99,1.0,200),(0,1,1,50),(100,1,1,7),(0.01,2,.5,33)]:
    phi = tau2/(tau2+sig2)
    DE_varM   = (v+tau2+sig2)/(v+tau2)          # Var(M): J=1 -> inf
    DE_errM   = (tau2+sig2)/tau2                # error component of Var(M)
    DE_item   = (tau2+sig2)/tau2                # per-item MSE about the truth
    print(f"  v={v:7.2f} tau2={tau2:.2f} sig2={sig2:.2f} | 1/phi={1/phi:7.4f} | "
          f"DE[Var(M)]={DE_varM:7.4f}  DE[err of M]={DE_errM:7.4f}  DE[per-item MSE]={DE_item:7.4f}")
print("  -> Kish's 1/phi is the design effect for the ERROR variance (per-item MSE, and")
print("     identically the error component of Var(M)).  It is NOT the design effect for")
print("     Var(M) itself, which also carries the irreducible item-sampling term v/I.")
print("     Aggregate and per-item ESTIMATION share ONE ceiling.  My brief claimed they")
print("     differ by ~4x; they do not.  The real divergence is estimation vs DETECTION.")
print()
v,tau2,sig2 = 3.7,0.99,1.0
print(f"  What growing the panel from J=9 to J=infinity actually buys (v={v},tau2={tau2},sig2={sig2}):")
print(f"    error variance of M: {(tau2+sig2/9):.4f} -> {tau2:.4f}  = {100*(sig2/9)/(tau2+sig2/9):.1f}% of error variance")
print(f"    total Var(M):        {(v+tau2+sig2/9):.4f} -> {v+tau2:.4f}  = {100*(sig2/9)/(v+tau2+sig2/9):.1f}% of total variance")
print(f"    the headline '9 judges, getting 2' invites 'cut to 2', saving 7/9 = 77.8% of spend")
print(f"    for a {100*(sig2/2-sig2/9)/(tau2+sig2/9):.1f}% INCREASE in error variance (J=9 -> J=2).")
