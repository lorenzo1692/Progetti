"""Fit of CASE_STRESS_SURROGATE (physics/case_stress_surrogate.m).

Target: maximum primary Pm and Pm+Pb over the case SCLs of the 2D FE
(validation/results/mech_reference_fe.csv, 19 runs). Base: membrane stress
intensity of the nose of the layered model, Pm_LC
(validation/results/case_surrogate_features.csv, run_case_surrogate_features.m).
Correction for the frame action the axisymmetric model misses:

    sigma = (a + b ln x) * Pm_LC,   x = h_WP / t_nose

Prints the coefficients, the fit error and the leave-one-design-out
(a design and its jacket-thickness variants left out together).
"""
import csv, collections, os
import numpy as np
from scipy.optimize import least_squares

here = os.path.dirname(os.path.abspath(__file__))
res = os.path.join(here, '..', 'results')
F = {r['name']: {k: (float(v) if k != 'name' else v) for k, v in r.items()}
     for r in csv.DictReader(open(os.path.join(res, 'case_surrogate_features.csv')))}
R = collections.defaultdict(dict)
for r in csv.DictReader(open(os.path.join(res, 'mech_reference_fe.csv'))):
    if r['kind'] == 'scl':
        R[r['name']][r['label']] = (float(r['Pm_P']), float(r['PmPb_P']))
names = list(F)
grp = [n.split('@')[0] for n in names]
x = np.array([F[n]['h_wp'] / F[n]['t_n'] for n in names])
base = np.array([F[n]['lcPm'] for n in names])
model = lambda c, xx: c[0] + c[1] * np.log(xx)
for j, lbl in ((0, 'Pm'), (1, 'Pm+Pb')):
    y = np.array([max(v[j] for v in R[n].values()) for n in names])
    fit = lambda m: least_squares(lambda c: model(c, x[m]) * base[m] / y[m] - 1, [1, 0.1]).x
    c = fit(np.ones(len(names), bool))
    r_all = model(c, x) * base / y
    loo = []
    for g in sorted(set(grp)):
        m = np.array([q != g for q in grp])
        loo += list(model(fit(m), x[~m]) * base[~m] / y[~m])
    loo = np.array(loo)
    print(f'{lbl:6s} a = {c[0]:.8f}, b = {c[1]:.8f} | fit rms {100*np.sqrt(np.mean((r_all-1)**2)):.1f} % | '
          f'LOO model/FE {loo.min():.3f} .. {loo.max():.3f}, rms {100*np.sqrt(np.mean((loo-1)**2)):.1f} %')
print(f'x range {x.min():.2f} .. {x.max():.2f}')
