import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
# receptor: 7FIH chain R protein atoms only
out=[]
for L in open('/mnt/user-data/uploads/7FIH.pdb'):
    if L[:4]=='ATOM' and L[21]=='R' and L[16] in (' ','A'):
        out.append(L)
open('/home/claude/dock/receptor.pdb','w').write(''.join(out)+'END\n')
print("receptor atoms:",len(out))

# ligand SDF (V2000) from coordinates + derived connectivity, with Kekule bond orders
names,bonds,E=pickle.load(open('/home/claude/super/bonds.pkl','rb'))
order={}
def setb(a,b,o): order[frozenset((a,b))]=o
for a in names:
    for b in bonds[a]: setb(a,b,1)
# Kekule assignments
for a,b in [('C8','C9'),('C10','C15'),('C11','N12'),('C13','N14'),
            ('C19','C20'),('C21','C22'),('C23','C24'),('C6','O7'),('C26','O27')]:
    setb(a,b,2)
lines=[]
lines.append("55Z_Org43553")
lines.append("  computed from 7FIH coordinates")
lines.append("")
bl=sorted({frozenset(k) for k in order})
lines.append(f"{len(names):3d}{len(bl):3d}  0  0  0  0  0  0  0  0999 V2000")
idx={n:i+1 for i,n in enumerate(names)}
for n in names:
    x,y,z=L55[n][0]
    lines.append(f"{x:10.4f}{y:10.4f}{z:10.4f} {L55[n][1]:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
for k in bl:
    a,b=tuple(k)
    lines.append(f"{idx[a]:3d}{idx[b]:3d}{order[k]:3d}  0  0  0  0")
lines.append("M  END"); lines.append("$$$$")
open('/home/claude/dock/xtal_lig.sdf','w').write("\n".join(lines)+"\n")
print("ligand sdf written:",len(names),"atoms",len(bl),"bonds")
# sanity: bond order sums per atom (valence check, heavy atoms only, ignoring H)
val={}
for k,o in order.items():
    a,b=tuple(k); val[a]=val.get(a,0)+o; val[b]=val.get(b,0)+o
maxv={'C':4,'N':3,'O':2,'S':6}
bad=[(n,val[n],E[names.index(n)]) for n in names if val[n]>maxv[E[names.index(n)]]]
print("valence violations:",bad if bad else "none")
