"""
gphr_lib.py -- shared primitives for the LHCGR/FSHR allosteric-pocket analysis.

INPUTS EXPECTED (paths passed by callers; defaults point at ../inputs/):
    7FIH.pdb   LHCGR(S277I) + CG + Gs + Org 43553 (ligand 55Z), chain R = receptor
    8I2G.pdb   FSHR + FSH + Gs + Cpd-21f  (ligand O6F), chain R = receptor

No third-party dependencies beyond numpy.
Hydrogens are ignored throughout. Altlocs other than ' ' and 'A' are skipped.
"""
import numpy as np

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

VDW = {'C':1.70,'N':1.55,'O':1.52,'S':1.80}
COV = {'C':0.77,'N':0.75,'O':0.73,'S':1.02}

_B62_ORDER = "ARNDCQEGHILKMFPSTWYV"
_B62_RAW = """
 4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0
-1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3
-2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3
-2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3
 0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1
-1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2
-1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2
 0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3
-2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3
-1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3
-1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1
-1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2
-1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1
-2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1
-1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2
 1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2
 0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0
-3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3
-2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1
 0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4
"""
_M = np.array([[int(v) for v in row.split()] for row in _B62_RAW.strip().split('\n')])
BLOSUM62 = {(a, b): int(_M[i, j])
            for i, a in enumerate(_B62_ORDER) for j, b in enumerate(_B62_ORDER)}


def parse_pdb(fn, chain):
    """Return (residues, order, hetero).

    residues : {(resnum, resname): {atomname: (xyz, element)}}   protein, given chain
    order    : [(resnum, resname)] in file order
    hetero   : {(resname, chain, resnum): {atomname: (xyz, element)}}
    """
    res, order, het = {}, [], {}
    for L in open(fn):
        rec = L[:6].strip()
        if rec not in ('ATOM', 'HETATM'):
            continue
        if L[16] not in (' ', 'A'):
            continue
        ch, rn, name = L[21], L[17:20].strip(), L[12:16].strip()
        num = int(L[22:26])
        el = (L[76:78].strip() or name[0]).upper()
        xyz = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
        if rec == 'ATOM' and ch == chain and rn in AA3:
            k = (num, rn)
            if k not in res:
                res[k] = {}
                order.append(k)
            res[k][name] = (xyz, el)
        elif rec == 'HETATM':
            het.setdefault((rn, ch, num), {})[name] = (xyz, el)
    return res, order, het


def sequence(order):
    return ''.join(AA3[k[1]] for k in order)


def needleman_wunsch(a, b, matrix=None, gap_open=-11, gap_extend=-1,
                     match=2, mismatch=-1):
    """Global alignment, affine gaps, three-state (Gotoh). Returns [(i, j)] pairs.

    matrix=BLOSUM62 uses substitution scores; matrix=None uses match/mismatch.
    """
    n, m = len(a), len(b)
    NEG = -1e9
    M = np.full((n + 1, m + 1), NEG)
    X = np.full((n + 1, m + 1), NEG)
    Y = np.full((n + 1, m + 1), NEG)
    M[0, 0] = 0.0
    for i in range(1, n + 1):
        X[i, 0] = gap_open + gap_extend * (i - 1)
    for j in range(1, m + 1):
        Y[0, j] = gap_open + gap_extend * (j - 1)
    pM = np.zeros((n + 1, m + 1), np.int8)
    pX = np.zeros((n + 1, m + 1), np.int8)
    pY = np.zeros((n + 1, m + 1), np.int8)
    for i in range(1, n + 1):
        ai = a[i - 1]
        for j in range(1, m + 1):
            bj = b[j - 1]
            if matrix is None:
                s = match if ai == bj else mismatch
            else:
                s = matrix.get((ai, bj), mismatch)
            c = (M[i - 1, j - 1], X[i - 1, j - 1], Y[i - 1, j - 1])
            k = int(np.argmax(c)); M[i, j] = c[k] + s; pM[i, j] = k
            c = (M[i - 1, j] + gap_open, X[i - 1, j] + gap_extend, Y[i - 1, j] + gap_open)
            k = int(np.argmax(c)); X[i, j] = c[k]; pX[i, j] = k
            c = (M[i, j - 1] + gap_open, X[i, j - 1] + gap_open, Y[i, j - 1] + gap_extend)
            k = int(np.argmax(c)); Y[i, j] = c[k]; pY[i, j] = k
    i, j = n, m
    state = int(np.argmax([M[n, m], X[n, m], Y[n, m]]))
    ptr = (pM, pX, pY)
    pairs = []
    while i > 0 or j > 0:
        if state == 0:
            pairs.append((i - 1, j - 1)); ns = ptr[0][i, j]; i -= 1; j -= 1; state = ns
        elif state == 1:
            ns = ptr[1][i, j]; i -= 1; state = ns
        else:
            ns = ptr[2][i, j]; j -= 1; state = ns
        if i == 0 and j == 0:
            break
        if i == 0:
            state = 2
        elif j == 0:
            state = 1
    return pairs[::-1]


def kabsch(P, Q):
    """Rotation R and translation t minimising |R P + t - Q|."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def heavy(atoms):
    return {a: v for a, v in atoms.items() if v[1] != 'H'}


def min_dist_to_ligand(atoms, lig_xyz):
    d = 1e9
    for a, (x, e) in atoms.items():
        if e == 'H':
            continue
        d = min(d, np.linalg.norm(lig_xyz - x, axis=1).min())
    return d


def shell(res, lig, cutoff):
    """Residues whose min heavy-atom distance to the ligand is < cutoff."""
    lx = np.array([v[0] for v in lig.values() if v[1] != 'H'])
    out = {}
    for k, ats in res.items():
        d = min_dist_to_ligand(ats, lx)
        if d < cutoff:
            out[k] = d
    return out


def derive_bonds(lig):
    """Covalent-radius connectivity. Returns {atom: [neighbours]}."""
    names = [n for n, v in lig.items() if v[1] != 'H']
    X = np.array([lig[n][0] for n in names])
    E = [lig[n][1] for n in names]
    bonds = {n: [] for n in names}
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            lim = COV.get(E[i], .8) + COV.get(E[j], .8) + 0.45
            if np.linalg.norm(X[i] - X[j]) < lim:
                bonds[names[i]].append(names[j])
                bonds[names[j]].append(names[i])
    return bonds, names, E


def neighbour_counts(lig, res, radius=5.0):
    """Protein heavy atoms within `radius` of each ligand heavy atom."""
    prot = np.array([v[0] for k, ats in res.items() for a, v in ats.items() if v[1] != 'H'])
    return {n: int((np.linalg.norm(prot - v[0], axis=1) < radius).sum())
            for n, v in lig.items() if v[1] != 'H'}


def write_pdb_ligand(path, lig, resname, chain='X', resnum=1):
    with open(path, 'w') as f:
        for i, (n, (x, e)) in enumerate(sorted(lig.items()), 1):
            f.write(f"HETATM{i:5d} {n:<4}{resname:>3} {chain}{resnum:4d}    "
                    f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          {e:>2}\n")
        f.write("END\n")
