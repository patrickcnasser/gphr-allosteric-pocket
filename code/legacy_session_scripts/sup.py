import numpy as np, sys
AA3={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def parse(fn, chain):
    res={}; order=[]
    het={}
    for L in open(fn):
        rec=L[:6].strip()
        if rec not in ('ATOM','HETATM'): continue
        ch=L[21]; rn=L[17:20].strip(); num=int(L[22:26]); name=L[12:16].strip()
        el=L[76:78].strip() or name[0]
        xyz=np.array([float(L[30:38]),float(L[38:46]),float(L[46:54])])
        alt=L[16]
        if alt not in (' ','A'): continue
        if rec=='ATOM' and ch==chain and rn in AA3:
            k=(num,rn)
            if k not in res: res[k]={}; order.append(k)
            res[k][name]=(xyz,el)
        if rec=='HETATM':
            het.setdefault((rn,ch,num),{})[name]=(xyz,el)
    return res,order,het

def seq(order): return ''.join(AA3[k[1]] for k in order)

def nw(a,b,match=2,mis=-1,gap=-8,ext=-1):
    n,m=len(a),len(b)
    NEG=-1e9
    M=np.full((n+1,m+1),NEG); X=np.full((n+1,m+1),NEG); Y=np.full((n+1,m+1),NEG)
    M[0,0]=0
    for i in range(1,n+1): X[i,0]=gap+ext*(i-1)
    for j in range(1,m+1): Y[0,j]=gap+ext*(j-1)
    ptrM=np.zeros((n+1,m+1),int); ptrX=np.zeros((n+1,m+1),int); ptrY=np.zeros((n+1,m+1),int)
    for i in range(1,n+1):
        for j in range(1,m+1):
            s=match if a[i-1]==b[j-1] else mis
            cand=[M[i-1,j-1],X[i-1,j-1],Y[i-1,j-1]]
            k=int(np.argmax(cand)); M[i,j]=cand[k]+s; ptrM[i,j]=k
            cand=[M[i-1,j]+gap, X[i-1,j]+ext, Y[i-1,j]+gap]
            k=int(np.argmax(cand)); X[i,j]=cand[k]; ptrX[i,j]=k
            cand=[M[i,j-1]+gap, X[i,j-1]+gap, Y[i,j-1]+ext]
            k=int(np.argmax(cand)); Y[i,j]=cand[k]; ptrY[i,j]=k
    i,j=n,m; state=int(np.argmax([M[n,m],X[n,m],Y[n,m]])); pairs=[]
    ptr=[ptrM,ptrX,ptrY]
    while i>0 or j>0:
        if state==0:
            pairs.append((i-1,j-1)); ns=ptr[0][i,j]; i-=1; j-=1; state=ns
        elif state==1:
            ns=ptr[1][i,j]; i-=1; state=ns
        else:
            ns=ptr[2][i,j]; j-=1; state=ns
        if i==0 and j==0: break
        if i==0: state=2
        elif j==0: state=1
    return pairs[::-1]

def kabsch(P,Q):
    pc=P.mean(0); qc=Q.mean(0)
    H=(P-pc).T@(Q-qc)
    U,S,Vt=np.linalg.svd(H)
    d=np.sign(np.linalg.det(Vt.T@U.T))
    D=np.diag([1,1,d])
    R=Vt.T@D@U.T
    return R,qc-R@pc

L_res,L_ord,L_het=parse('/mnt/user-data/uploads/7FIH.pdb','R')
F_res,F_ord,F_het=parse('/mnt/user-data/uploads/8I2G.pdb','R')
lig55Z={k:v for k,v in L_het.items() if k[0]=='55Z'}
ligO6F={k:v for k,v in F_het.items() if k[0]=='O6F'}
L55=list(lig55Z.values())[0]; O6F=list(ligO6F.values())[0]
L55={n:v for n,v in L55.items() if v[1]!='H'}
O6F={n:v for n,v in O6F.items() if v[1]!='H'}
print("55Z heavy atoms:",len(L55)," O6F heavy atoms:",len(O6F))
print("LHCGR chain R residues:",len(L_ord),"range",L_ord[0][0],"-",L_ord[-1][0])
print("FSHR  chain R residues:",len(F_ord),"range",F_ord[0][0],"-",F_ord[-1][0])

sa=seq(L_ord); sb=seq(F_ord)
pairs=nw(sa,sb)
ident=sum(1 for i,j in pairs if sa[i]==sb[j])
print(f"aligned pairs {len(pairs)}  identity {100*ident/len(pairs):.1f}%")
amap={L_ord[i][0]:F_ord[j][0] for i,j in pairs}
for q in (612,515,585,589,350,451,595,349,519,531,593,604):
    if q in amap:
        fi=[k for k in F_ord if k[0]==amap[q]][0]
        li=[k for k in L_ord if k[0]==q][0]
        print(f"  LHCGR {li[1]}{q} -> FSHR {fi[1]}{amap[q]}")
np.save('/home/claude/super/pairs.npy',np.array(pairs))
import pickle
pickle.dump((L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap),open('/home/claude/super/data.pkl','wb'))
