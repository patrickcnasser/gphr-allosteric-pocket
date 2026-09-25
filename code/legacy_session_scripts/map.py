import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
def kabsch(P,Q):
    pc=P.mean(0); qc=Q.mean(0); H=(P-pc).T@(Q-qc)
    U,S,Vt=np.linalg.svd(H); d=np.sign(np.linalg.det(Vt.T@U.T))
    R=Vt.T@np.diag([1,1,d])@U.T; return R,qc-R@pc

ligc=np.mean([v[0] for v in L55.values()],axis=0)
# TM-bundle selection: aligned pairs whose LHCGR CA is within 25 A of the ligand centroid AND resnum>=340 (TMD)
P=[];Q=[];lab=[]
for i,j in pairs:
    lk=L_ord[i]; fk=F_ord[j]
    if 'CA' not in L_res[lk] or 'CA' not in F_res[fk]: continue
    ca=L_res[lk]['CA'][0]
    if lk[0]<340: continue
    if np.linalg.norm(ca-ligc)>25: continue
    P.append(ca); Q.append(F_res[fk]['CA'][0]); lab.append((lk,fk))
P=np.array(P);Q=np.array(Q)
print("TM-core CA pairs selected:",len(P))
# iterative outlier pruning
keep=np.ones(len(P),bool)
for it in range(8):
    R,t=kabsch(P[keep],Q[keep])
    d=np.linalg.norm((P@R.T+t)-Q,axis=1)
    rms=np.sqrt((d[keep]**2).mean())
    cut=max(2.0, d[keep].mean()+2*d[keep].std())
    nk=d<cut
    print(f"  iter {it}: n={keep.sum()} RMSD={rms:.2f} A  cutoff={cut:.2f}")
    if (nk==keep).all(): break
    keep=nk
R,t=kabsch(P[keep],Q[keep])
d=np.linalg.norm((P@R.T+t)-Q,axis=1)
print(f"FINAL: {keep.sum()} CA pairs, RMSD {np.sqrt((d[keep]**2).mean()):.2f} A")

# map 55Z into FSHR frame
M55={n:(v[0]@R.T+t, v[1]) for n,v in L55.items()}
# sanity: how well does mapped 55Z overlap the real FSHR ligand O6F?
o=np.array([v[0] for v in O6F.values()])
m=np.array([v[0] for v in M55.values()])
print(f"centroid offset mapped-55Z vs O6F: {np.linalg.norm(m.mean(0)-o.mean(0)):.2f} A")
# fraction of mapped 55Z atoms within 2 A of any O6F atom
dd=np.linalg.norm(m[:,None,:]-o[None,:,:],axis=2)
print(f"mapped 55Z atoms within 2.0 A of an O6F atom: {(dd.min(1)<2.0).sum()}/{len(m)}")

def contacts(lig, resd, cut=4.5):
    out={}
    for k,ats in resd.items():
        best=1e9; batom=None; lat=None
        for an,(xyz,el) in ats.items():
            if el=='H': continue
            for ln,(lxyz,lel) in lig.items():
                dist=np.linalg.norm(xyz-lxyz)
                if dist<best: best=dist; batom=an; lat=ln
        if best<cut: out[k]=(best,batom,lat)
    return out

cL=contacts(L55,L_res)
cF=contacts(M55,F_res)
print(f"\nLHCGR native contacts <4.5 A: {len(cL)}")
print(f"FSHR  mapped contacts <4.5 A: {len(cF)}")

# clashes in FSHR
print("\n--- CLASHES: FSHR heavy atoms < 3.2 A of mapped 55Z ---")
rows=[]
for k,ats in F_res.items():
    for an,(xyz,el) in ats.items():
        if el=='H': continue
        for ln,(lxyz,lel) in M55.items():
            dist=np.linalg.norm(xyz-lxyz)
            if dist<3.2: rows.append((dist,k,an,ln))
rows.sort()
seen=set()
for dist,k,an,ln in rows:
    if k in seen: continue
    seen.add(k)
    # is this position divergent?
    inv={v:kk for kk,v in amap.items()}
    lres=inv.get(k[0])
    lname=[x[1] for x in L_ord if x[0]==lres]
    lname=lname[0] if lname else '?'
    div='DIVERGENT' if lname!=k[1] else 'conserved'
    print(f"  {dist:.2f} A  FSHR {k[1]}{k[0]} ({an}) vs 55Z {ln}   [LHCGR {lname}{lres}] {div}")

# per-residue comparison at the equivalent positions
print("\n--- CONTACT COMPARISON at aligned positions ---")
inv={v:kk for kk,v in amap.items()}
allpos=set()
for k in cL: allpos.add(('L',k))
print(f"{'LHCGR res':<12}{'dist':>7}   {'FSHR res':<12}{'dist':>7}   delta")
for k in sorted(cL, key=lambda x: cL[x][0]):
    fnum=amap.get(k[0])
    fk=[x for x in F_ord if x[0]==fnum]
    if not fk: 
        print(f"  {k[1]}{k[0]:<8}{cL[k][0]:>7.2f}   (unaligned)")
        continue
    fk=fk[0]
    fd=cF.get(fk,(None,))[0]
    ds=f"{fd:>7.2f}" if fd else "   >4.5"
    dv='*' if fk[1]!=k[1] else ' '
    delta=f"{fd-cL[k][0]:+.2f}" if fd else ""
    print(f"{dv} {k[1]}{k[0]:<9}{cL[k][0]:>7.2f}   {fk[1]}{fk[0]:<9}{ds}   {delta}")
pickle.dump((M55,R,t,cL,cF,amap),open('/home/claude/super/mapped.pkl','wb'))
