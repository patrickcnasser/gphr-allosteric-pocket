import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
M55,R,t,cL,cF,amap=pickle.load(open('/home/claude/super/mapped.pkl','rb'))
# forward was x_F = R x_L + t ; inverse x_L = R^T (x_F - t)
MO6={n:(R.T@(v[0]-t), v[1]) for n,v in O6F.items()}
print("=== control: round-trip the forward map ===")
back={n:(R.T@(v[0]-t)) for n,v in M55.items()}
err=max(np.linalg.norm(back[n]-L55[n][0]) for n in L55)
print(f"  max round-trip error on 55Z: {err:.4f} A (should be ~0)")

# where does mapped O6F sit relative to 55Z?
o=np.array([v[0] for v in MO6.values()]); x=np.array([v[0] for v in L55.values()])
print(f"  centroid offset mapped-O6F vs 55Z in LHCGR frame: {np.linalg.norm(o.mean(0)-x.mean(0)):.2f} A")

print("\n=== CLASH SCAN: LHCGR heavy atoms vs mapped Cpd-21f ===")
prot=[(k,an,v[0]) for k,ats in L_res.items() for an,v in ats.items() if v[1]!='H']
px=np.array([p[2] for p in prot])
rows=[]
for ln,(xx,e) in MO6.items():
    d=np.linalg.norm(px-xx,axis=1); i=int(d.argmin())
    rows.append((d[i],prot[i][0],prot[i][1],ln))
rows.sort()
nsev=sum(1 for r in rows if r[0]<2.5); ntight=sum(1 for r in rows if r[0]<3.0)
for d,k,an,ln in rows[:12]:
    flag='  <-- SEVERE' if d<2.5 else ('  <-- tight' if d<3.0 else '')
    print(f"  {d:.2f} A  LHCGR {k[1]}{k[0]} ({an})  vs O6F {ln}{flag}")
print(f"\n  atoms <3.0 A: {ntight}/{len(rows)}   atoms <2.5 A: {nsev}/{len(rows)}")

print("\n=== THE TEST: does Tyr612 clash with Cpd-21f? ===")
lk=[k for k in L_res if k[0]==612][0]
best=(9e9,None,None)
for an,(xx,e) in L_res[lk].items():
    if e=='H': continue
    for ln,(yy,e2) in MO6.items():
        d=np.linalg.norm(xx-yy)
        if d<best[0]: best=(d,an,ln)
print(f"  closest LHCGR Tyr612 - Cpd-21f contact: {best[0]:.2f} A  ({best[1]} ... {best[2]})")
# compare with what His615 does in the real FSHR structure
fk=[k for k in F_res if k[0]==615][0]
b2=(9e9,None,None)
for an,(xx,e) in F_res[fk].items():
    if e=='H': continue
    for ln,(yy,e2) in O6F.items():
        d=np.linalg.norm(xx-yy)
        if d<b2[0]: b2=(d,an,ln)
print(f"  native FSHR His615 - Cpd-21f contact:    {b2[0]:.2f} A  ({b2[1]} ... {b2[2]})")
print(f"  => Tyr612 is {best[0]-b2[0]:+.2f} A relative to the native histidine contact")

# per-Tyr612-atom detail
print("\n  Tyr612 atom-by-atom vs nearest Cpd-21f atom:")
for an in ['CB','CG','CD1','CD2','CE1','CE2','CZ','OH']:
    if an in L_res[lk]:
        d=min(np.linalg.norm(L_res[lk][an][0]-v[0]) for v in MO6.values())
        print(f"    {an:<4} {d:.2f} A")

print("\n=== which LHCGR residues are the real problem, if any ===")
seen={}
for d,k,an,ln in rows:
    if d<4.0 and k not in seen: seen[k]=(d,an,ln)
inv={v:kk for kk,v in amap.items()}
for k,(d,an,ln) in sorted(seen.items(), key=lambda x:x[1][0])[:10]:
    fn=amap.get(k[0]); fr=[y for y in F_ord if y[0]==fn]
    dv='DIVERGENT' if (fr and fr[0][1]!=k[1]) else 'conserved'
    fname=f"{fr[0][1]}{fn}" if fr else '?'
    print(f"  {d:.2f} A  {k[1]}{k[0]} ({an}) vs O6F {ln}   [FSHR {fname}] {dv}")
