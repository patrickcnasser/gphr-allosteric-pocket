import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
M55,R,t,cL,cF,amap=pickle.load(open('/home/claude/super/mapped.pkl','rb'))
rad={'C':0.77,'N':0.75,'O':0.73,'S':1.02}
nm=list(O6F); X=np.array([O6F[n][0] for n in nm]); E=[O6F[n][1] for n in nm]
b={n:[] for n in nm}
for i in range(len(nm)):
    for j in range(i+1,len(nm)):
        if np.linalg.norm(X[i]-X[j])<rad[E[i]]+rad[E[j]]+0.45:
            b[nm[i]].append(nm[j]); b[nm[j]].append(nm[i])
# find tert-butyl: quaternary C with 3 terminal C
tb=[n for n in nm if E[nm.index(n)]=='C' and len(b[n])==4 and sum(1 for x in b[n] if len(b[x])==1)==3]
print("O6F quaternary carbons:",tb)
for q in tb:
    grp=[q]+[x for x in b[q] if len(b[x])==1]
    print("  tert-butyl group:",grp)
    hits=[]
    for k,at in F_res.items():
        d=min(np.linalg.norm(v[0]-O6F[a][0]) for an,v in at.items() if v[1]!='H' for a in grp)
        if d<4.5: hits.append(f"{k[1]}{k[0]}({d:.1f})")
    print("   FSHR residues <4.5 A of O6F tert-butyl:"," ".join(hits))
# and where is 55Z tert-butyl mapped relative to O6F tert-butyl
tb55=['C1','C2','C3','C4']
c1=np.mean([M55[a][0] for a in tb55],axis=0)
for q in tb:
    grp=[q]+[x for x in b[q] if len(b[x])==1]
    c2=np.mean([O6F[a][0] for a in grp],axis=0)
    print(f"  centroid offset 55Z tert-butyl (mapped) vs O6F tert-butyl: {np.linalg.norm(c1-c2):.2f} A")
# His615 rotamer check: chi1
fk=[k for k in F_res if k[0]==615][0]; fa=F_res[fk]
print("\nHis615 atoms present:", sorted(a for a in fa if fa[a][1]!='H'))
lk=[k for k in L_res if k[0]==612][0]
print("Tyr612 atoms present:", sorted(a for a in L_res[lk] if L_res[lk][a][1]!='H'))
# distance from His615 ring atoms to mapped 55Z tert-butyl
for an in ['CG','ND1','CD2','CE1','NE2']:
    if an in fa:
        d=min(np.linalg.norm(fa[an][0]-M55[a][0]) for a in tb55)
        print(f"  His615 {an} -> nearest 55Z tert-butyl atom: {d:.2f} A")
for an in ['CG','CD1','CD2','CE1','CE2','CZ','OH']:
    if an in L_res[lk]:
        d=min(np.linalg.norm(L_res[lk][an][0]-L55[a][0]) for a in tb55)
        print(f"  Tyr612 {an} -> nearest 55Z tert-butyl atom: {d:.2f} A")
