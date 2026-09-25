import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
M55,R,t,cL,cF,amap=pickle.load(open('/home/claude/super/mapped.pkl','rb'))
# derive bonds by covalent radii
rad={'C':0.77,'N':0.75,'O':0.73,'S':1.02}
names=list(L55.keys()); X=np.array([L55[n][0] for n in names]); E=[L55[n][1] for n in names]
bonds={n:[] for n in names}
for i in range(len(names)):
    for j in range(i+1,len(names)):
        d=np.linalg.norm(X[i]-X[j]); lim=rad[E[i]]+rad[E[j]]+0.45
        if d<lim: bonds[names[i]].append(names[j]); bonds[names[j]].append(names[i])
print("formula check: C%d N%d O%d S%d"%(E.count('C'),E.count('N'),E.count('O'),E.count('S')))
print("expected Org 43553 C24 N6 O3 S2 (heavy=35)\n")
for n in names: print(f"  {n:<5}{E[names.index(n)]}  deg {len(bonds[n])}  -> {','.join(sorted(bonds[n]))}")
pickle.dump((names,bonds,E),open('/home/claude/super/bonds.pkl','wb'))
