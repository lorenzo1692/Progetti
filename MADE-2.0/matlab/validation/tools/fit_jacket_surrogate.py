"""Fit and leave-one-design-out check of JACKET_STRESS_SURROGATE.

Data: ../results/jacket_surrogate_calibration.csv, one row per layer, from
validation/jacket_surrogate_fe_data.m (2D FE WP_MECH_SURROGATE, primary load
case). Model (see physics/jacket_stress_surrogate.m):

    stress_k = S_z + a*sigma_nom_k + b*sigma_acc_k

sigma_acc_k = q_k/p_rs*sigma_nom_k, q_k = accumulated n*Iop*B up to the middle
of layer k per unit first-layer width. Prints the coefficients for Pm and
Pm+Pb, the leave-one-design-out errors, and the same errors for the analytic
formula of the scan (SCF table x SCF_transition_provisional x dcr_WP_rad).
"""
import csv, os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
R = list(csv.DictReader(open(os.path.join(here, '..', 'results', 'jacket_surrogate_calibration.csv'))))
f = lambda k: np.array([float(r[k]) for r in R])
name = np.array([r['name'] for r in R]); k = f('k').astype(int)
W = f('W'); Iop = f('Iop'); nt = f('n_turns'); Bl = f('Bl'); T = f('T'); param = f('param')
sig = f('sigma_nom'); Sz = f('S_z'); dcr = f('dcr_WP_rad'); PmPb = f('PmPb_P'); Pm = f('Pm_P')
rst = f('r_steel'); dj = f('dcr_jckt')
p_rs = sig/(rst*dj)
sacc = np.zeros_like(sig)
for n in set(name):
    m = np.where(name == n)[0]; m = m[np.argsort(k[m])]
    F = nt[m]*Iop[m]*Bl[m]
    q = (np.cumsum(F) - 0.5*F)/W[m][0]
    sacc[m] = q/p_rs[m]*sig[m]
designs = sorted(set(name))
X = np.column_stack([sig, sacc])

def report(label, ref, pred):
    e = pred/ref - 1
    emax = [pred[name == d].max()/ref[name == d].max() - 1 for d in designs]
    print('%-40s layer rms %5.1f%%  max %5.1f%% | design max %+5.1f .. %+5.1f %%'
          % (label, 100*np.sqrt((e**2).mean()), 100*abs(e).max(), 100*min(emax), 100*max(emax)))

for label, ref in [('Pm', Pm), ('Pm+Pb', PmPb)]:
    y = ref - Sz
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = np.zeros_like(y)
    for d in designs:
        te = name == d
        cd, *_ = np.linalg.lstsq(X[~te], y[~te], rcond=None)
        pred[te] = X[te] @ cd
    print('%s: a = %.3f, b = %.3f  (%d layers, %d designs)' % (label, c[0], c[1], len(y), len(designs)))
    report('  surrogate, leave-one-design-out', ref, Sz + pred)
    xxx = [1.087,1.136,1.190,1.250,1.316,1.389,1.471,1.563,1.667,1.786,1.923,2.083,2.273,2.500,2.778]
    yyy = [1.01,1.03,1.06,1.10,1.16,1.22,1.29,1.36,1.43,1.50,1.57,1.64,1.71,1.78,1.85]
    scf = np.polyval(np.polyfit(xxx, yyy, 5), param)
    report('  analytic formula of the scan', ref, Sz + dcr*sig*scf*np.where(T > 0, 3.15, 1))
