import numpy as np, pickle
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
M55,R,t,cL,cF,amap=pickle.load(open('/home/claude/super/bonds.pkl','rb')) if False else pickle.load(open('/home/claude/super/mapped.pkl','rb'))
TB=['C1','C2','C3','C4']   # tert-butyl: C2 quaternary, C1/C3/C4 methyls
vdw={'C':1.70,'N':1.55,'O':1.52,'S':1.80}

def cavity(lig, resd, anchor_atoms, probe=1.4, box=8.0, step=0.5):
    """Grid points that are (a) empty of protein, (b) empty of ligand, (c) contiguous-ish near the anchor."""
    c=np.mean([lig[a][0] for a in anchor_atoms],axis=0)
    prot=[(v[0],vdw.get(v[1],1.7)) for k,ats in resd.items() for an,v in ats.items() if v[1]!='H']
    px=np.array([p[0] for p in prot]); pr=np.array([p[1] for p in prot])
    m=np.linalg.norm(px-c,axis=1)<box+12; px=px[m]; pr=pr[m]
    lx=np.array([v[0] for v in lig.values()]); lr=np.array([vdw.get(v[1],1.7) for v in lig.values()])
    g=np.arange(-box,box+1e-9,step)
    G=np.stack(np.meshgrid(g,g,g,indexing='ij'),-1).reshape(-1,3)+c
    # keep only points within box radius
    G=G[np.linalg.norm(G-c,axis=1)<=box]
    keep=[]
    for i in range(0,len(G),4000):
        ch=G[i:i+4000]
        dp=np.linalg.norm(ch[:,None,:]-px[None,:,:],axis=2)-pr[None,:]
        dl=np.linalg.norm(ch[:,None,:]-lx[None,:,:],axis=2)-lr[None,:]
        ok=(dp>probe).all(1)&(dl>probe).all(1)&((dp<6.0).sum(1)>25)
        keep.append(ch[ok])
    cav=np.concatenate(keep) if keep else np.zeros((0,3))
    return c,cav

print("=== empty, buried space adjacent to the tert-butyl ===")
for tag,lig,resd in [('LHCGR (7FIH)',L55,L_res),('FSHR  (mapped)',M55,F_res)]:
    c,cav=cavity(lig,resd,TB)
    if len(cav)==0:
        print(f"{tag}: no buried empty volume found adjacent to tert-butyl"); continue
    d=np.linalg.norm(cav-c,axis=1)
    print(f"{tag}: {len(cav)} probe-accessible grid points, volume ~{len(cav)*0.4**3:.0f} A^3, max reach {d.max():.1f} A from tert-butyl centroid")
    # what lines the far end
    far=cav[d>d.max()-1.5]
    lin={}
    for k,ats in resd.items():
        m=min(np.linalg.norm(v[0]-far,axis=1).min() for an,v in ats.items() if v[1]!='H')
        if m<5.0: lin[f"{k[1]}{k[0]}"]=m
    print("   lined by:", ", ".join(f"{k}({v:.1f})" for k,v in sorted(lin.items(),key=lambda x:x[1])))

print("\n=== H-bond partner geometry from the tert-butyl ===")
c_tb_L=np.mean([L55[a][0] for a in TB],axis=0)
c_tb_F=np.mean([M55[a][0] for a in TB],axis=0)
lk=[k for k in L_res if k[0]==612][0]; fk=[k for k in F_res if k[0]==615][0]
yOH=L_res[lk]['OH'][0]
print(f"LHCGR Tyr612 OH  -> tert-butyl centroid: {np.linalg.norm(yOH-c_tb_L):.2f} A")
for an in ['ND1','NE2']:
    print(f"FSHR  His615 {an} -> tert-butyl centroid: {np.linalg.norm(F_res[fk][an][0]-c_tb_F):.2f} A")
# quaternary carbon as the growth origin
print(f"\nfrom the quaternary carbon C2:")
print(f"  LHCGR Tyr612 OH : {np.linalg.norm(yOH-L55['C2'][0]):.2f} A")
for an in ['ND1','NE2','CE1']:
    print(f"  FSHR His615 {an}: {np.linalg.norm(F_res[fk][an][0]-M55['C2'][0]):.2f} A")
# amide N as origin (Moore's N-methyl position)
print(f"\nfrom the amide nitrogen N5 (Moore compound 7 position):")
print(f"  LHCGR Tyr612 OH : {np.linalg.norm(yOH-L55['N5'][0]):.2f} A")
for an in ['ND1','NE2']:
    print(f"  FSHR His615 {an}: {np.linalg.norm(F_res[fk][an][0]-M55['N5'][0]):.2f} A")
# other divergent partners near tert-butyl
print("\n=== the other two divergent residues at this fragment ===")
for lnum,fnum in [(515,518),(531,534)]:
    a=[k for k in L_res if k[0]==lnum][0]; b=[k for k in F_res if k[0]==fnum][0]
    dl=min(np.linalg.norm(v[0]-L55[x][0]) for an,v in L_res[a].items() if v[1]!='H' for x in TB)
    df=min(np.linalg.norm(v[0]-M55[x][0]) for an,v in F_res[b].items() if v[1]!='H' for x in TB)
    print(f"  LHCGR {a[1]}{lnum} {dl:.2f} A  |  FSHR {b[1]}{fnum} {df:.2f} A   delta {df-dl:+.2f}")
