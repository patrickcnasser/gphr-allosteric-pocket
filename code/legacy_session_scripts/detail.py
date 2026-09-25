import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
M55,R,t,cL,cF,amap=pickle.load(open('/home/claude/super/mapped.pkl','rb'))
def get(resd,num):
    for k in resd:
        if k[0]==num: return k,resd[k]
    return None,None

print("=== the two lost contacts, atom by atom ===")
for lnum,fnum,tag in [(612,615,'position 7.42'),(595,598,'tail lysine')]:
    lk,la=get(L_res,lnum); fk,fa=get(F_res,fnum)
    # LHCGR: closest side-chain atom to native ligand
    best=(1e9,None,None)
    for an,(x,e) in la.items():
        if e=='H' or an in ('N','CA','C','O'): continue
        for ln,(lx,le) in L55.items():
            d=np.linalg.norm(x-lx)
            if d<best[0]: best=(d,an,ln)
    bestF=(1e9,None,None)
    for an,(x,e) in fa.items():
        if e=='H' or an in ('N','CA','C','O'): continue
        for ln,(lx,le) in M55.items():
            d=np.linalg.norm(x-lx)
            if d<bestF[0]: bestF=(d,an,ln)
    print(f"{tag}: LHCGR {lk[1]}{lnum} {best[1]}-{best[2]} {best[0]:.2f} A  |  FSHR {fk[1]}{fnum} {bestF[1]}-{bestF[2]} {bestF[0]:.2f} A")
    # CA-CA after superposition to see if it's sidechain identity or backbone shift
    ca_l=la['CA'][0]@R.T+t; ca_f=fa['CA'][0]
    print(f"    CA-CA after superposition: {np.linalg.norm(ca_l-ca_f):.2f} A")

print("\n=== what does FSHR His615 do instead? distance to real FSHR ligand O6F ===")
fk,fa=get(F_res,615)
b=(1e9,None,None)
for an,(x,e) in fa.items():
    if e=='H': continue
    for ln,(lx,le) in O6F.items():
        d=np.linalg.norm(x-lx)
        if d<b[0]: b=(d,an,ln)
print(f"  His615 {b[1]} - O6F {b[2]}: {b[0]:.2f} A")

print("\n=== burial: protein heavy-atom neighbours within 5 A per ligand atom ===")
def nb(lig,resd):
    prot=np.array([v[0] for k,ats in resd.items() for an,v in ats.items() if v[1]!='H'])
    out={}
    for ln,(x,e) in lig.items():
        out[ln]=int((np.linalg.norm(prot-x,axis=1)<5.0).sum())
    return out
nL=nb(L55,L_res); nF=nb(M55,F_res)
tot_l=sum(nL.values()); tot_f=sum(nF.values())
print(f"  total neighbour count  LHCGR {tot_l}   FSHR(mapped) {tot_f}   ratio {tot_f/tot_l:.2f}")
worst=sorted(((nF[k]-nL[k],k) for k in nL))
print("  biggest per-atom losses:", [(k,nL[k],nF[k]) for d,k in worst[:5]])
print("  biggest per-atom gains :", [(k,nL[k],nF[k]) for d,k in worst[-5:]])

print("\n=== minimum FSHR distance for every 55Z atom (clash scan) ===")
prot=[(k,an,v[0]) for k,ats in F_res.items() for an,v in ats.items() if v[1]!='H']
px=np.array([p[2] for p in prot])
bad=0
for ln,(x,e) in sorted(M55.items()):
    d=np.linalg.norm(px-x,axis=1); i=int(d.argmin())
    flag=''
    if d[i]<3.0: flag='  <-- tight'; bad+=1
    print(f"  {ln:<5} {d[i]:.2f} A to {prot[i][0][1]}{prot[i][0][0]} {prot[i][1]}{flag}")
print(f"\natoms under 3.0 A: {bad}")
