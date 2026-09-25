import numpy as np, pickle, sys
sys.path.insert(0,'/home/claude/super'); sys.path.insert(0,'/home/claude/species')
from sup import nw, AA3
from seqs import SEQ
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
hum=''.join(AA3[k[1]] for k in L_ord)

# CONTROL: does fetched human P22888 contain the structure sequence?
p=nw(hum,SEQ['human_P22888'])
ident=sum(1 for i,j in p if hum[i]==SEQ['human_P22888'][j])
off=set(L_ord[i][0]-(j+1) for i,j in p)
print(f"CONTROL human P22888 vs 7FIH chain R: {100*ident/len(p):.1f}% identity over {len(p)} pairs")
print(f"  numbering offsets present: {sorted(off)}  (0 = UniProt numbering == structure numbering)")
mis=[(L_ord[i][0],L_ord[i][1],SEQ['human_P22888'][j]) for i,j in p if hum[i]!=SEQ['human_P22888'][j]]
print(f"  mismatches: {mis if mis else 'none'}")

lig=np.array([v[0] for v in L55.values()])
shell=[]
for k,ats in L_res.items():
    d=min(np.linalg.norm(v[0]-lig,axis=1).min() for an,v in ats.items() if v[1]!='H')
    if d<8.0: shell.append((k,d))
shell.sort(key=lambda x:x[1])
contact=[s for s in shell if s[1]<4.5]
print(f"\npocket shells: {len(shell)} within 8 A, {len(contact)} within 4.5 A")

print("\n=== POSITION 515 ACROSS SPECIES (aligned, not eyeballed) ===")
res={}
for name,s in SEQ.items():
    pr=nw(hum,s)
    m={L_ord[i][0]:s[j] for i,j in pr}
    res[name]=m
    print(f"  {name:<20} 515-equivalent = {m.get(515,'?')}")
print(f"  FSHR 518 = L   TSHR 570 = L  (from structure/literature)")

print("\n=== FULL 8 A POCKET: divergences from human, per species ===")
for name,m in res.items():
    if name=='human_P22888': continue
    div=[(k[0],k[1],m.get(k[0],'?'),d) for k,d in shell if m.get(k[0]) and m[k[0]]!=AA3[k[1]]]
    print(f"\n  {name}: {len(div)} of {len(shell)} divergent")
    for n,h,o,d in sorted(div,key=lambda x:x[3]):
        mark='  <== CONTACT SHELL' if d<4.5 else ''
        print(f"     {h}{n} -> {o}   (min dist {d:.2f} A){mark}")
