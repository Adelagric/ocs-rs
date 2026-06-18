"""
Prototype: a "support-first" exact solver for genomic Optimum Contribution
Selection, exploring whether the ultra-sparsity of the OCS solution can break
the O(n^3) wall of the conic interior-point method (Clarabel).

This is a research prototype (numpy/scipy), not production code. It backs the
note in SUPPORT_FIRST.md. Run:  python3 research/support_first.py

Idea: the OCS optimum activates only a tiny support S (a handful of candidates).
For a fixed support, eliminating the equality 1'c=1 turns the problem into
"maximise a linear form over an ellipsoid" -> closed form (two linear solves).
The whole cost is then identifying S, done by column generation: add by reduced
cost / diversify by least relatedness, all via matrix-vector products G c that
never form G (G c = ridge*c + Z(Z' c)/s). Cost ~ |S| products of O(n*m), versus
~17 dense KKT factorisations of O(n^3) for the IPM.

Verified here: exactness vs scipy/Clarabel to 1e-8..1e-12; |S| stays bounded
(does not grow with n) even under a structured pedigree, correlated EBV and a
tight kinship bound.
"""

import numpy as np
from scipy.linalg import null_space

RIDGE = 1e-5

def sim_pop(n,m,seed,structured,n_founders=24,n_gen=8):
    rng=np.random.default_rng(seed); p=rng.uniform(0.05,0.5,m)
    if not structured:
        M=rng.binomial(2,p,size=(n,m)).astype(float)
    else:
        cur=rng.binomial(2,p,size=(n_founders,m)).astype(float)   # fondateurs
        for g in range(n_gen):                                     # ségrégation mendélienne
            n_off = n if g==n_gen-1 else max(n_founders, n//3)
            sires=rng.integers(0,len(cur),n_off); dams=rng.integers(0,len(cur),n_off)
            cur=(rng.binomial(1,cur[sires]/2)+rng.binomial(1,cur[dams]/2)).astype(float)
        M=cur
    Z=M-2*p; s=2*np.sum(p*(1-p))
    return Z,s,p

def make_ebv(Z,seed,correlated):
    rng=np.random.default_rng(seed+999)
    if not correlated: return rng.standard_normal(Z.shape[0])
    g=Z@rng.standard_normal(Z.shape[1])                           # valeur génétique additive
    return (g-g.mean())/g.std()

def matvec(Z,s,c): return RIDGE*c+(Z@(Z.T@c))/s
def GSS(Z,s,S): ZS=Z[S]; return ZS@ZS.T/s+RIDGE*np.eye(len(S))
def closed_form(Z,s,b,k,S):
    GS=GSS(Z,s,S); bS=b[S]; nS=len(S)
    if nS==1: return (np.array([1.]),bS[0],0.) if GS[0,0]<=k else None
    c0=np.ones(nS)/nS; N=null_space(np.ones((1,nS)))
    Gt=N.T@GS@N; bt=N.T@bS; gt=N.T@(GS@c0); q0=c0@GS@c0
    gi=np.linalg.solve(Gt,gt); bi=np.linalg.solve(Gt,bt); rho2=gt@gi-(q0-k)
    if rho2<=0: return None
    cS=c0+N@(-gi+np.sqrt(rho2)*bi/np.sqrt(bt@bi))
    (mu,lam),*_=np.linalg.lstsq(np.column_stack([np.ones(nS),2*(GS@cS)]),bS,rcond=None)
    return cS,mu,lam
def support_first(Z,s,b,k,tol=1e-7,max_iter=4000):
    n=len(b); nprod=0; S=[int(np.argmax(b))]; c_cur=np.zeros(n); c_cur[S[0]]=1.
    for it in range(max_iter):
        sol=closed_form(Z,s,b,k,S)
        if sol is None:
            Gc=matvec(Z,s,c_cur); nprod+=1
            for j in np.argsort(Gc):
                if j not in S: S.append(int(j)); break
            continue
        cS,mu,lam=sol
        if cS.min()<-tol: S=[S[i] for i in range(len(S)) if cS[i]>tol]; continue
        c=np.zeros(n); c[S]=cS; c_cur=c; Gc=matvec(Z,s,c); nprod+=1
        r=b-mu-2*lam*Gc; r[S]=-np.inf; j=int(np.argmax(r))
        if r[j]<=tol: return c,sorted(S),it+1,nprod
        S.append(j)
    return c_cur,sorted(S),max_iter,nprod

def krange(Z,s):
    n=Z.shape[0]; d=RIDGE+(Z*Z).sum(1)/s
    one=np.ones(n); Gone=matvec(Z,s,one); kmin=(one@Gone)/n**2   # uniforme
    return kmin, d.mean()                                          # (serré, ~lâche)

# exactitude vs scipy sur pop structurée+corrélée (n=400)
from scipy.optimize import minimize
Z,s,p=sim_pop(400,5000,1,True); b=make_ebv(Z,1,True); G=Z@Z.T/s+RIDGE*np.eye(400)
kmin,kmax=krange(Z,s)
print("=== exactitude (pédigree + EBV corrélé, n=400) ===")
for frac in [0.2,0.5,0.9]:
    k=kmin+frac*(kmax-kmin); c,S,it,npd=support_first(Z,s,b,k)
    cons=[{"type":"eq","fun":lambda c:c.sum()-1},{"type":"ineq","fun":lambda c:k-c@G@c}]
    r=minimize(lambda c:-(b@c),np.ones(400)/400,bounds=[(0,None)]*400,constraints=cons,
               method="SLSQP",options={"ftol":1e-11,"maxiter":3000})
    print(f"  k={k:.4f}  |S|={len(S):>3}  gain_sf={b@c:.5f}  gain_scipy={-r.fun:.5f}  Δ={abs(b@c+r.fun):.1e}")

# stress |S|(k) : iid vs pédigree, EBV iid vs corrélé, n=2000
print("\n=== |S| en fonction de k (n=2000, m=20000) ===")
print(f"  {'régime':<28} " + " ".join(f"k={f:>4.2f}" for f in [0.05,0.15,0.30,0.50,0.80]))
for tag,struct,corr in [("iid geno + iid EBV",False,False),
                        ("iid geno + EBV corrélé",False,True),
                        ("pédigree + iid EBV",True,False),
                        ("pédigree + EBV corrélé",True,True)]:
    Z,s,p=sim_pop(2000,20000,3,struct); b=make_ebv(Z,3,corr)
    kmin,kmax=krange(Z,s)
    row=[]
    for frac in [0.05,0.15,0.30,0.50,0.80]:
        k=kmin+frac*(kmax-kmin); c,S,it,npd=support_first(Z,s,b,k)
        row.append(f"{len(S):>5}")
    print(f"  {tag:<28} " + " ".join(row) + f"   (kmin={kmin:.3f} kmax={kmax:.3f})")

print("\n=== |S|(n) a k SERRE fixe, pedigree + EBV correle (le pire regime) ===")
print(f"  {'frac_k':>7} " + " ".join(f"n={n:>5}" for n in [500,2000,5000,10000]))
for frac in [0.05, 0.15, 0.30]:
    cells=[]
    for n in [500,2000,5000,10000]:
        Z,s,p=sim_pop(n,20000,3,True); b=make_ebv(Z,3,True)
        kmin,kmax=krange(Z,s); k=kmin+frac*(kmax-kmin)
        c,S,it,npd=support_first(Z,s,b,k)
        cells.append(f"|S|={len(S):>4}({npd:>3}p)")
    print(f"  {frac:>7.2f} " + " ".join(cells))
print("  (p = produits matrice-vecteur O(nm) = cout dominant)")
