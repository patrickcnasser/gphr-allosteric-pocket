import numpy as np, pickle, sys
sys.path.insert(0,'/home/claude/super')
from sup import nw, AA3
L_res,L_ord,F_res,F_ord,L55,O6F,pairs,amap=pickle.load(open('/home/claude/super/data.pkl','rb'))
rat="".join("""MGRRVPALRQLLVLAVLLLKPSQLQSRELSGSRCPEPCDCAPDGALRCPGPRAGLARLSL
TYLPVKVIPSQAFRGLNEVVKIEISQSDSLERIEANAFDNLLNLSELLIQNTKNLLYIEP
GAFTNLPRLKYLSICNTGIRTLPDVTKISSSEFNFILEICDNLHITTIPGNAFQGMNNES
VTLKLYGNGFEEVQSHAFNGTTLISLELKENIYLEKMHSGAFQGATGPSILDISSTKLQA
LPSHGLESIQTLIALSSYSLKTLPSKEKFTSLLVATLTYPSHCCAFRNLPKKEQNFSFSI
FENFSKQCESTVRKADNETLYSAIFEENELSGWDYDYGFCSPKTLQCAPEPDAFNPCEDI
MGYAFLRVLIWLINILAIFGNLTVLFVLLTSRYKLTVPRFLMCNLSFADFCMGLYLLLIA
SVDSQTKGQYYNHAIDWQTGSGCGAAGFFTVFASELSVYTLTVITLERWHTITYAVQLDQ
KLRLRHAIPIMLGGWLFSTLIATMPLVGISNYMKVSICLPMDVESTLSQVYILSILILNV
VAFVVICACYIRIYFAVQNPELTAPNKDTKIAKKMAILIFTDFTCMAPISFFAISAAFKV
PLITVTNSKILLVLFYPVNSCANPFLYAIFTKAFQRDFLLLLSRFGCCKRRAELYRRKEF
SAYTSNCKNGFPGASKPSQATLKLSTVHCQQPIPPRALTH""".split())
print("rat P16235 length:",len(rat),"(UniProt states 700)")
hum=''.join(AA3[k[1]] for k in L_ord)
pairs2=nw(hum,rat)
ident=sum(1 for i,j in pairs2 if hum[i]==rat[j])
print(f"human-structure vs rat: {len(pairs2)} aligned pairs, {100*ident/len(pairs2):.1f}% identity")
ratmap={L_ord[i][0]:(j+1,rat[j],hum[i]) for i,j in pairs2}

# pocket shells
lig=np.array([v[0] for v in L55.values()])
def shell(cut):
    s=[]
    for k,ats in L_res.items():
        d=min(np.linalg.norm(v[0]-lig,axis=1).min() for an,v in ats.items() if v[1]!='H')
        if d<cut: s.append((k,d))
    return sorted(s,key=lambda x:x[1])
for cut,label in [(8.0,'8 A shell'),(4.5,'4.5 A contact shell')]:
    sh=shell(cut)
    div=[]
    miss=[]
    for k,d in sh:
        if k[0] not in ratmap: miss.append(k); continue
        rn,ra,ha=ratmap[k[0]]
        if ra!=ha: div.append((k,d,rn,ra))
    print(f"\n{label}: {len(sh)} residues, {len(div)} divergent rat vs human, {len(miss)} unaligned")
    for k,d,rn,ra in div:
        print(f"   human {k[1]}{k[0]} (min dist {d:.2f} A)  ->  rat {ra}{rn}")
    if not div: print("   ALL IDENTICAL")
# explicit check of the named contact residues
print("\nnamed residues:")
for n in (350,451,515,519,528,531,585,589,593,595,604,612,349,597,601,592,599,447,450,535,588,517,516):
    if n in ratmap:
        rn,ra,ha=ratmap[n]
        print(f"   {ha}{n} -> rat {ra}{rn} {'DIVERGENT' if ra!=ha else ''}")
