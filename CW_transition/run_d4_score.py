"""Score the pre-registered D4 realisation-level discordance gate from the banks. CPU only.

Run from the repo root:  python3 CW_transition/run_d4_score.py
Banks: reports/ignite_bank.npz (IGNITE Stage-1), reports/ig_fnull*.npz (IGNITE-2 fresh nulls),
reports/ig_sloop*.npz + ig_sloopX*.npz (IGNITE-2 soft loop), reports/ig2_floors.npz.
Writes reports/d4_score.npz.
"""
import numpy as np, glob, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from criterion_v2_1_d4 import Levers, S_det_of, d4_gate, PREREG, CRIT2

ALPHA = 0.10          # T-CAL false-alarm rate on true-signal realisations. FIXED IN ADVANCE.
FORMS = ['oracle', 'obs']
lev = Levers()
B = np.load('reports/ignite_bank.npz', allow_pickle=True)
F = np.load('reports/ig2_floors.npz', allow_pickle=True)

CELLS = [(-12.75, 30, 'lit', float(F['fN_0']), float(F['fN_sd_0'])),
         (-13.25, 30, 'vlbi', float(F['fN_1']), float(F['fN_sd_1']))]


def bank_pop(kind, h, T, tier, floor):
    idx = np.where((B['kind'] == kind) & (B['h'] == h) & (B['T_label'] == T)
                   & (B['tier'] == tier) & (B['tol'] == 0.0))[0]
    out = []
    for i in idx:
        rec = {f: S_det_of(lev, B['mapk'][i] - B['n_true_grid'][i], B['mapk'][i], B['qmax'][i],
                           B['dlnL_det'][i], B['lnK'][i], tier, int(B['T_label'][i]),
                           int(B['geo_id'][i]), floor, f) for f in FORMS}
        rec['tag'] = '%s#%d' % (kind, i)
        rec['dk'] = B['mapk'][i] - B['n_true_grid'][i]
        rec['names'] = B['names']
        out.append(rec)
    return out


def npz_pop(pattern, floor, it=None):
    out = []
    for fn in sorted(glob.glob(pattern)):
        d = np.load(fn, allow_pickle=True)
        tier = str(d['tier']); T = int(d['T_label']); geo = int(d['geo_id'])
        if it is None:
            mapk, qm, dlnL = d['mapk'], d['qmax'], d['dlnL_det']
        else:
            k = it if it >= 0 else int(d['n_iter']) - 1
            mapk, qm, dlnL = d['iter_mapk'][k], d['iter_qmax'][k], d['iter_dlnL'][k]
        dk = mapk - d['n_true_grid']
        rec = {f: S_det_of(lev, dk, mapk, qm, dlnL, d['lnK'], tier, T, geo, floor, f)
               for f in FORMS}
        rec['tag'] = os.path.basename(fn)
        rec['dk'] = dk
        rec['names'] = d['names']
        out.append(rec)
    return out


def detecting(pop, form):
    return [r for r in pop if r[form]['n_det'] > 0]


def flagrate(pop, form, thr):
    d = detecting(pop, form)
    if not d:
        return np.nan, 0, 0
    f = sum(1 for r in d if np.isfinite(r[form]['S']) and r[form]['S'] > thr)
    return f / len(d), f, len(d)


print(PREREG)
print('=' * 100)
print('VALUE GATE (the D4 machinery must reproduce IGNITE-2\'s banked statistic first):')
gate_ok = d4_gate()
print('  gate: %s' % ('ALL PASS' if gate_ok else 'FAIL -- SCORING ABORTED'))
if not gate_ok:
    sys.exit(1)

results = {}
for h, T, tier, floor, fsd in CELLS:
    hs = '%d' % round(abs(h) * 100)
    cell = '(%s, %d, %s)' % (h, T, tier)
    print('\n' + '=' * 100)
    print('CELL %s   fresh v2 Gumbel floor fN = %.2f +- %.2f nat (alpha=0.05, N=150)' % (cell, floor, fsd))

    true_cal = bank_pop('sig', h, T, tier, floor)                                  # calibration
    true_test = npz_pop('reports/ig_sloop_h%s_T%d_%s_*.npz' % (hs, T, tier), floor, it=0)
    true_test_f = npz_pop('reports/ig_sloop_h%s_T%d_%s_*.npz' % (hs, T, tier), floor, it=-1)
    wrong = (bank_pop('nullA', h, T, tier, floor) + bank_pop('nullL', h, T, tier, floor)
             + npz_pop('reports/ig_fnullA_h%s_T%d_%s_*.npz' % (hs, T, tier), floor)
             + npz_pop('reports/ig_fnullL_h%s_T%d_%s_*.npz' % (hs, T, tier), floor))
    scram = npz_pop('reports/ig_sloopX_h%s_T%d_%s_*.npz' % (hs, T, tier), floor, it=0)
    scram_f = npz_pop('reports/ig_sloopX_h%s_T%d_%s_*.npz' % (hs, T, tier), floor, it=-1)
    noise = bank_pop('nullN', h, T, tier, floor) + npz_pop('reports/ig_fnullN_h%s_T%d_%s_*.npz' % (hs, T, tier), floor)

    for form in FORMS:
        cal = np.array([r[form]['S'] for r in detecting(true_cal, form)], float)
        cal = cal[np.isfinite(cal)]
        t_cal = float(np.quantile(cal, 1 - ALPHA)) if len(cal) else np.nan
        print('\n  --- form = D4-%s ---' % form.upper())
        print('  T-CAL calibration set: %d detecting true-signal reals (IGNITE Stage-1); '
              'S_det median=%.3g  q90=%.3g' % (len(cal), np.median(cal), t_cal))
        for tname, thr in [('T-CHI2', CRIT2), ('T-CAL', t_cal)]:
            r_i, fi, ni = flagrate(wrong, form, thr)
            r_s, fs, ns = flagrate(scram, form, thr)
            allw_f, allw_n = fi + fs, ni + ns
            r_ii, f2, n2 = flagrate(true_test, form, thr)
            r_iif, f2f, n2f = flagrate(true_test_f, form, thr)
            r_cal, fc, nc = flagrate(true_cal, form, thr)
            r_iv, f4, n4 = flagrate(noise, form, thr)
            n_vac = sum(1 for r in noise if r[form]['n_det'] == 0)
            print('   [%s thr=%.4g]' % (tname, thr))
            print('     (i)   wrong-counterpart flagged : %3d/%-3d = %5.1f%%   (nullA/L+fresh %d/%d, scrambled-loop %d/%d)   %s'
                  % (allw_f, allw_n, 100 * allw_f / max(allw_n, 1), fi, ni, fs, ns,
                     'PASS' if allw_n and allw_f / allw_n >= 0.95 else 'FAIL'))
            print('     (ii)  true-signal flagged  OOS  : %3d/%-3d = %5.1f%%   (sloop it0)   %s'
                  % (f2, n2, 100 * r_ii if n2 else float("nan"),
                     'PASS' if n2 and r_ii <= 0.10 else 'FAIL'))
            print('           true-signal flagged  OOS  : %3d/%-3d = %5.1f%%   (sloop fixed point)'
                  % (f2f, n2f, 100 * r_iif if n2f else float('nan')))
            print('           true-signal flagged in-sample (calibration set) : %3d/%-3d = %5.1f%%'
                  % (fc, nc, 100 * r_cal if nc else float('nan')))
            print('     (iv)  pure-noise: %d/%d vacuous (no detection) -> NO FLAG by construction; '
                  'of the %d that DO detect, flagged %d = %.1f%%'
                  % (n_vac, len(noise), n4, f4, 100 * r_iv if n4 else float('nan')))
            results['%s|%s|%s' % (cell, form, tname)] = np.array(
                [thr, allw_f, allw_n, f2, n2, f2f, n2f, fc, nc, f4, n4, n_vac], float)

    # ---- (iii) the small-|dk| keepers: scrambled-loop realisations that KEEP a certification
    print('\n  --- (iii) THE SMALL-|dk| KEEPERS (the D1 hole; per-pulsar form missed these) ---')
    for r in scram_f:
        st = r['oracle']
        if st['n_cert'] == 0:
            continue
        ci = np.where(st['cert'])[0]
        for form in FORMS:
            cal = np.array([q[form]['S'] for q in detecting(true_cal, form)], float)
            cal = cal[np.isfinite(cal)]
            t_cal = float(np.quantile(cal, 1 - ALPHA))
            S = r[form]['S']
            print('    %-42s cert=%-11s dk=%-5d [%s] S_det=%10.4g  chi2bar:%-6s calbar(%.4g):%s'
                  % (r['tag'][:42], str(r['names'][ci[0]]), int(r['dk'][ci[0]]), form.upper()[:3],
                     S, 'FLAG' if S > CRIT2 else 'miss', t_cal, 'FLAG' if S > t_cal else 'MISS'))

np.savez('reports/d4_score.npz', alpha=ALPHA, crit2=CRIT2, prereg=PREREG,
         keys=np.array(list(results.keys())),
         vals=np.array(list(results.values())) if results else np.zeros((0, 12)),
         cols=np.array(['thr', 'wrong_flag', 'wrong_n', 'true_oos_flag', 'true_oos_n',
                        'true_fix_flag', 'true_fix_n', 'true_cal_flag', 'true_cal_n',
                        'noise_flag', 'noise_n', 'noise_vacuous']))
print('\nbanked -> reports/d4_score.npz')
