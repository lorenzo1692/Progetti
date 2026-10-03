#!/usr/bin/env python3
"""Peak field per turn and per layer from an ANSYS 2D EM listing export.

Usage:  python3 em_layer_peaks.py <TFBM_em_*.txt> [B_peak_layers comma list]

The export (e.g. TFBM_em_plane233.txt, conductor-only PLANE233 model) must
contain, in this order: NLIST (node coordinates), ELIST (elements with MAT,
TYPE, REAL: one real set per turn) and PRNSOL,B (BX BY BZ BSUM). Layers are
the groups of turns with the same centroid height y (tolerance 5 mm). If the
per-layer peaks of the MADE scan (column B_peak_layers of the results file)
are given, they are compared layer by layer.
"""
import re
import sys
import numpy as np


def parse(path):
    nodes, elems, B = {}, [], {}
    mode = None
    with open(path, errors='replace') as f:
        for ln in f:
            if 'LIST ALL SELECTED NODES' in ln:
                mode = 'n'; continue
            if 'LIST ALL SELECTED ELEMENTS' in ln:
                mode = 'e'; continue
            if 'PRINT B ' in ln or 'PRINT B\t' in ln:
                mode = 'b'; continue
            if 'PRINT ' in ln or 'ELEMENT TABLE LISTING' in ln:
                mode = None; continue
            s = ln.split()
            if not s or not re.match(r'^\d+$', s[0]):
                continue
            try:
                if mode == 'n' and len(s) >= 4:
                    nodes[int(s[0])] = (float(s[1]), float(s[2]))
                elif mode == 'e' and len(s) >= 8:
                    elems.append([int(x) for x in s])
                elif mode == 'b' and len(s) == 5:
                    B[int(s[0])] = float(s[4])
            except ValueError:
                pass
    return nodes, elems, B


def main():
    nodes, elems, B = parse(sys.argv[1])
    turn = {}
    for e in elems:
        for n in e[6:]:
            if n in B:
                turn.setdefault(e[3], set()).add(n)
    T = []
    for r, ns in sorted(turn.items()):
        xy = np.array([nodes[n] for n in ns])
        T.append((r, max(B[n] for n in ns), xy[:, 0].mean(), xy[:, 1].mean()))
    T = np.array(T)
    o = np.argsort(-T[:, 3])
    layers = []
    for i in o:
        if layers and abs(T[layers[-1][0], 3] - T[i, 3]) < 5e-3:
            layers[-1].append(i)
        else:
            layers.append([i])
    ref = [float(v) for v in sys.argv[2].split(',')] if len(sys.argv) > 2 else None
    print('%d nodes, %d elements, %d turns, peak BSUM on conductor %.3f T' % (len(nodes), len(elems), len(T), max(B.values())))
    print('layer  y [m]   turns  ANSYS peak [T]' + ('   MADE [T]   diff [T]' if ref else ''))
    for k, L in enumerate(layers):
        pk = T[L, 1].max()
        line = 'L%-4d %.4f %5d   %8.3f' % (k + 1, T[L, 3].mean(), len(L), pk)
        if ref and k < len(ref):
            line += '   %8.3f   %+7.3f' % (ref[k], ref[k] - pk)
        print(line)


if __name__ == '__main__':
    main()
