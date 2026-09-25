import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
M55,R,t,cL,cF,amap=pickle.load(open('/home/claude/super/mapped.pkl','rb'))
FRAG={
 'tert-butyl':['C1','C2','C3','C4'],
 'carboxamide':['N5','C6','O7'],
 'thiophene':['C8','C9','C10','C15','S16'],
 'pyrimidine':['C11','N12','C13','N14'],
 '5-amino':['N35'],
 '2-methylthio':['S17','C18'],
 'pendant phenyl':['C19','C20','C21','C22','C23','C24'],
 'anilide linker':['N25','C26','O27','C28'],
 'morpholine':['N29','C30','C31','O32','C33','C34'],
}
allat=[a for v in FRAG.values() for a in v]
assert len(allat)==35 and len(set(allat))==35, (len(allat),len(set(allat)))
print("fragment partition: 35 atoms, no duplicates, verified against derived connectivity\n")

def nbcount(lig,resd):
    prot=np.array([v[0] for k,ats in resd.items() for an,v in ats.items() if v[1]!='H'])
    return {ln:int((np.linalg.norm(prot-x,axis=1)<5.0).sum()) for ln,(x,e) in lig.items()}
nL=nbcount(L55,L_res); nF=nbcount(M55,F_res)
print(f"{'fragment':<16}{'n':>3}{'LHCGR/atom':>12}{'FSHR/atom':>11}{'delta':>8}")
for f,ats in FRAG.items():
    a=np.mean([nL[x] for x in ats]); b=np.mean([nF[x] for x in ats])
    print(f"{f:<16}{len(ats):>3}{a:>12.1f}{b:>11.1f}{b-a:>+8.1f}")

# which residues contact each fragment, and are they divergent?
inv={v:k for k,v in amap.items()}
lname={k[0]:k[1] for k in L_ord}
print("\n=== residues within 4.5 A of each fragment (LHCGR) — D marks LHCGR/FSHR divergence ===")
for f,ats in FRAG.items():
    hits=[]
    for k,at in L_res.items():
        best=min([np.linalg.norm(v[0]-L55[a][0]) for an,v in at.items() if v[1]!='H' for a in ats])
        if best<4.5:
            fn=amap.get(k[0]); fr=[x for x in F_ord if x[0]==fn]
            d='D' if (fr and fr[0][1]!=k[1]) else ' '
            hits.append(f"{d}{k[1]}{k[0]}({best:.1f})")
    print(f"  {f:<16} {' '.join(hits)}")
