"""
rmsd_lib.py -- symmetry-corrected, order-independent heavy-atom RMSD for SDF poses.

WHY THIS EXISTS
gnina/openbabel do NOT preserve the input atom order when they write poses. Comparing
pose coordinates to the reference row-by-row therefore compares the wrong atoms and
produces nonsense (for this ligand, ~8 A where the true value is ~1 A).

This module establishes the correspondence from the molecular graph instead:
  1. read atoms and bonds from each SDF record;
  2. find a graph isomorphism between pose and reference over heavy atoms;
  3. enumerate the reference molecule's automorphisms (phenyl 2-fold, tert-butyl
     3-fold, morpholine 2-fold, ...) and report the minimum RMSD over all of them.

No third-party dependencies beyond numpy.
"""
import numpy as np
from itertools import permutations


def read_sdf(fn):
    """[(coords, elements, bonds)] over all records; bonds are 1-indexed (a, b, order)."""
    out, lines, i = [], open(fn).read().split('\n'), 0
    while i < len(lines):
        if i + 3 >= len(lines) or not lines[i + 3].strip():
            break
        try:
            na = int(lines[i + 3][:3]); nb = int(lines[i + 3][3:6])
        except ValueError:
            break
        xyz, el = [], []
        for j in range(i + 4, i + 4 + na):
            xyz.append([float(lines[j][0:10]), float(lines[j][10:20]), float(lines[j][20:30])])
            el.append(lines[j][31:34].strip())
        bonds = []
        for j in range(i + 4 + na, i + 4 + na + nb):
            bonds.append((int(lines[j][:3]), int(lines[j][3:6]), int(lines[j][6:9])))
        out.append((np.array(xyz), el, bonds))
        k = i
        while k < len(lines) and lines[k].strip() != '$$$$':
            k += 1
        i = k + 1
    return out


def heavy_graph(coords, el, bonds):
    """Drop hydrogens, reindex 0..n-1. Returns (coords, elements, adjacency sets)."""
    keep = [i for i, e in enumerate(el) if e != 'H']
    remap = {o: n for n, o in enumerate(keep)}
    adj = {i: set() for i in range(len(keep))}
    for a, b, _o in bonds:
        a -= 1; b -= 1
        if a in remap and b in remap:
            adj[remap[a]].add(remap[b]); adj[remap[b]].add(remap[a])
    return coords[keep], [el[i] for i in keep], adj


def _signatures(el, adj, rounds=4):
    """Morgan-style refinement: element + sorted neighbour signatures, iterated."""
    sig = {i: el[i] for i in adj}
    for _ in range(rounds):
        sig = {i: sig[i] + "(" + ",".join(sorted(sig[j] for j in adj[i])) + ")" for i in adj}
        # compress to keep strings short
        vocab = {s: str(k) for k, s in enumerate(sorted(set(sig.values())))}
        sig = {i: vocab[s] for i, s in sig.items()}
    return sig


def _match(adj_a, el_a, adj_b, el_b, limit=200000):
    """Yield isomorphisms a->b as dicts. Backtracking with signature pruning."""
    sa, sb = _signatures(el_a, adj_a), _signatures(el_b, adj_b)
    cand = {i: [j for j in adj_b if sb[j] == sa[i] and len(adj_b[j]) == len(adj_a[i])]
            for i in adj_a}
    if any(not v for v in cand.values()):
        return
    order = sorted(adj_a, key=lambda i: len(cand[i]))
    n = len(order)
    mapping, used = {}, set()
    steps = [0]

    def bt(k):
        if steps[0] > limit:
            return
        if k == n:
            yield dict(mapping)
            return
        i = order[k]
        for j in cand[i]:
            steps[0] += 1
            if j in used:
                continue
            ok = True
            for nb in adj_a[i]:
                if nb in mapping and mapping[nb] not in adj_b[j]:
                    ok = False; break
            if ok:
                for nb in adj_b[j]:
                    if nb in used and _inv(mapping, nb) not in adj_a[i]:
                        ok = False; break
            if not ok:
                continue
            mapping[i] = j; used.add(j)
            yield from bt(k + 1)
            del mapping[i]; used.discard(j)

    def _inv(m, v):
        for kk, vv in m.items():
            if vv == v:
                return kk
        return None

    yield from bt(0)


def symmetry_rmsd(pose_sdf_record, ref_sdf_record, max_autos=2000):
    """Minimum heavy-atom RMSD over graph isomorphism + reference automorphisms.

    Returns (rmsd, n_automorphisms, matched_bool).
    """
    Pc, Pe, Pb = pose_sdf_record
    Rc, Re, Rb = ref_sdf_record
    P, pe, padj = heavy_graph(Pc, Pe, Pb)
    R, re_, radj = heavy_graph(Rc, Re, Rb)
    if len(P) != len(R):
        return None, 0, False
    iso = next(_match(padj, pe, radj, re_), None)
    if iso is None:
        return None, 0, False
    # P reindexed into reference order
    Pm = np.zeros_like(R)
    for i, j in iso.items():
        Pm[j] = P[i]
    autos = []
    for a in _match(radj, re_, radj, re_):
        autos.append(a)
        if len(autos) >= max_autos:
            break
    best = None
    for a in autos:
        idx = [a[i] for i in range(len(R))]
        d = np.sqrt(((Pm[idx] - R) ** 2).sum(1).mean())
        if best is None or d < best:
            best = float(d)
    return best, len(autos), True


def centroid_offset(pose_rec, ref_rec):
    P, _pe, _ = heavy_graph(*pose_rec)
    R, _re, _ = heavy_graph(*ref_rec)
    return float(np.linalg.norm(P.mean(0) - R.mean(0)))


def atom_overlap(pose_rec, ref_rec, cut=2.0):
    P, _pe, _ = heavy_graph(*pose_rec)
    R, _re, _ = heavy_graph(*ref_rec)
    return float((np.linalg.norm(P[:, None, :] - R[None, :, :], axis=2).min(1) < cut).mean())


def mapped_coords(pose_rec, ref_rec):
    """Pose coordinates reordered into the reference atom order (None if no match)."""
    P, pe, padj = heavy_graph(*pose_rec)
    R, re_, radj = heavy_graph(*ref_rec)
    iso = next(_match(padj, pe, radj, re_), None)
    if iso is None:
        return None
    out = np.zeros_like(R)
    for i, j in iso.items():
        out[j] = P[i]
    return out
