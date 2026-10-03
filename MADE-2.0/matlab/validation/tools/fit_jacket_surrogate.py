"""Fit and leave-one-design-out check of JACKET_STRESS_SURROGATE.

Data: ../results/jacket_surrogate_calibration.csv, one row per layer, from
validation/jacket_surrogate_fe_data.m (2D FE WP_MECH_SURROGATE, primary load
case) and validation/run_jacket_jt_calibration.m (jacket-thickness variants,
names '<base>@jt+...'). Model (see physics/jacket_stress_surrogate.m):

    stress_k = S_z + (JT_k/JT_ref)^beta * (a*sigma_nom_k + b*sigma_acc_k [+ c*sigma_acc_k*s_k])

sigma_acc_k = q_k/p_rs*sigma_nom_k, q_k = accumulated n*Iop*B up to the middle
of layer k per unit first-layer width; s_k = (n_k - n_{k+1})/n_k the width
step under layer k (turns lost to the next layer). beta = 0 is the original
model (JT enters only through sigma_nom ~ 1/JT); beta > 0 weakens the JT
dependence. For every candidate form the script prints the leave-one-design-
out errors (all variants of a base design are left out together) and the
margin the scan needs so that the FE design maximum is covered:
margin = max over designs of FE_max / surrogate_max.

Usage: python3 fit_jacket_surrogate.py
"""
import csv
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
# Jacket-thickness variants whose WP no longer fits the sector (smallest
# lateral gap WP - case flank below 5 mm; the scan requires toroidal_gap =
# 15 mm and would have removed turns): the side wall becomes a thin
# ligament and the FE stress rises with JT instead of falling. Kept in the
# CSV, excluded from the fit (gap from jacket_jt_variant.m).
EXCLUDE = {'s1_1@jt+1.0': 'gap 2.0 mm', 's2_4@jt+1.0': 'gap -0.6 mm'}
R = [r for r in csv.DictReader(open(os.path.join(here, '..', 'results', 'jacket_surrogate_calibration.csv')))
     if r['name'] not in EXCLUDE]
f = lambda key: np.array([float(r[key]) for r in R])
name = np.array([r['name'] for r in R])
base = np.array([n.split('@')[0] for n in name])
k = f('k').astype(int)
W = f('W'); Iop = f('Iop'); nt = f('n_turns'); Bl = f('Bl'); T = f('T'); param = f('param')
sig = f('sigma_nom'); Sz = f('S_z'); dcr = f('dcr_WP_rad'); PmPb = f('PmPb_P'); Pm = f('Pm_P')
rst = f('r_steel'); dj = f('dcr_jckt'); JT_Ch = f('JT_Ch'); Ch_Cw = f('Ch_Cw')
p_rs = sig/(rst*dj)
sacc = np.zeros_like(sig); step = np.zeros_like(sig); JT = np.zeros_like(sig)
for n in set(name):
    m = np.where(name == n)[0]; m = m[np.argsort(k[m])]
    F = nt[m]*Iop[m]*Bl[m]
    q = (np.cumsum(F) - 0.5*F)/W[m][0]
    sacc[m] = q/p_rs[m]*sig[m]
    nn = np.r_[nt[m][1:], nt[m][-1]]
    step[m] = np.maximum(0, nt[m] - nn)/nt[m]
    Cw = W[m][0]/nt[m][0]
    JT[m] = JT_Ch[m]*Ch_Cw[m]*Cw
designs = sorted(set(base))
names = sorted(set(name))
JT_ref = 3e-3


def loo(X, y):
    pred = np.zeros_like(y)
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    for d in designs:
        te = base == d
        cd, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        pred[te] = X[te] @ cd
    return c, pred


def report(label, ref, pred):
    e = pred/ref - 1
    r = np.array([ref[name == n].max()/pred[name == n].max() for n in names])
    print('  %-34s layer rms %5.1f%% | design max %+6.1f .. %+6.1f %% | margin needed %.3f'
          % (label, 100*np.sqrt((e**2).mean()), 100*(1/r.max() - 1), 100*(1/r.min() - 1), r.max()))
    return r


def forms(beta):
    g = (JT/JT_ref)**beta
    return {
        'beta=%.2f' % beta: np.column_stack([g*sig, g*sacc]),
        'beta=%.2f + step' % beta: np.column_stack([g*sig, g*sacc, g*sacc*step]),
    }


print('%d layers, %d runs, %d base designs: %s' % (len(R), len(names), len(designs), ', '.join(designs)))
best = {}
for label, ref in [('Pm', Pm), ('Pm+Pb', PmPb)]:
    print('\n== %s ==' % label)
    y = ref - Sz
    res = []
    for beta in np.arange(0, 1.01, 0.1):
        for fl, X in forms(beta).items():
            c, pred = loo(X, y)
            e = (Sz + pred)/ref - 1
            rr = np.array([ref[name == n].max()/(Sz + pred)[name == n].max() for n in names])
            res.append((np.sqrt((e**2).mean()), fl, c, rr.max(), beta, X))
    res.sort(key=lambda t: t[0])
    for rms, fl, c, mg, beta, X in res[:6]:
        print('  %-18s rms %5.1f%%  margin %.3f  coef %s' % (fl, 100*rms, mg, np.round(c, 4)))
    c0, p0 = loo(forms(0)['beta=0.00'], y)
    print(' original form (beta = 0), coef %s:' % np.round(c0, 4)); report('original', ref, Sz + p0)
    rms, fl, c, mg, beta, X = res[0]
    _, pb = loo(X, y)
    print(' best form %s, coef %s:' % (fl, np.round(c, 4))); report('best', ref, Sz + pb)
    best[label] = (fl, beta, c)
    xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778]
    yyy = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85]
    scf = np.polyval(np.polyfit(xxx, yyy, 5), param)
    report('analytic formula of the scan', ref, Sz + dcr*sig*scf*np.where(T > 0, 3.15, 1))

# JT response on the variants: FE vs surrogate, design maximum relative to the base design
print('\n== JT response (design max Pm+Pb, ratio to the base design) ==')
fl, beta, c = best['Pm+Pb']
X = forms(beta)[fl]
_, pb = loo(X, PmPb - Sz)
_, p0 = loo(forms(0)['beta=0.00'], PmPb - Sz)
for d in designs:
    group = [n for n in names if n.split('@')[0] == d]
    if len(group) < 2:
        continue
    b0 = (name == d)
    for n in group:
        m = name == n
        print('  %-20s JT %.2f mm  FE %5.0f (%.3f)  best %5.0f (%.3f)  original %5.0f (%.3f)'
              % (n, 1e3*JT[m].mean(), PmPb[m].max()/1e6, PmPb[m].max()/PmPb[b0].max(),
                 (Sz + pb)[m].max()/1e6, (Sz + pb)[m].max()/(Sz + pb)[b0].max(),
                 (Sz + p0)[m].max()/1e6, (Sz + p0)[m].max()/(Sz + p0)[b0].max()))
