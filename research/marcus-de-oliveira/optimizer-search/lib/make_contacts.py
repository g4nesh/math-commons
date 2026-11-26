"""Builds contact_n4.json. Usage: python -I -B make_contacts.py <path to the audit moc4lib.py>"""
import sys, json, importlib.util
sys.path.insert(0, '.')
import numpy as np
from moc_hp import hp_margin
from moc_core import perms
from calibrate import load_config
spec = importlib.util.spec_from_file_location("m4", sys.argv[1])
m4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(m4)
P4 = [tuple(int(x) for x in p) for p in perms(4)]
assert [tuple(p) for p in m4.PERMS] == P4, "permutation order mismatch"
j = json.load(open('calibration/calib_n4c.json'))
num = []
for r in j['best']:
    if r['dist_genperm'] > 0.1 and r['dist_blockdiag'] > 0.1 and min(r['sep']) > 0.1:
        a, b, U = load_config(r)
        h = hp_margin(a, b, U, dps=60)
        num.append(dict(source="toolkit calibration run n4c (seed 105, restart %d)" % r['restart'],
                        kind="vertex contact" if r['n_active'] == 1 else "edge contact",
                        float_m=r['m'], hp_m_60digits=h['m_str'], active_perm_indices=r['active'],
                        active_perms=[list(P4[k]) for k in r['active']],
                        dist_genperm=r['dist_genperm'], dist_blockdiag=r['dist_blockdiag'], sep=r['sep'],
                        a=r['a'], b=r['b'], U_re=r['U_re'], U_im=r['U_im']))
out = dict(
  description="n = 4 contact data for the MOC search. Permutations are 0-based tuples in lexicographic order (index k = moc_core.perms(4)[k]); sigma[i] is the image of i.",
  explicit_configurations_from_audit=[],
  explicit_configurations_note="The Zheng n = 4 audit folder contains no explicit (a, b, U) configurations at which d touches the boundary of H; it holds only the combinatorial contact sets below.",
  combinatorial_contact_sets=dict(
    source="D. Zheng, The Marcus-de Oliveira conjecture for n = 4 (Oct 2026), Appendix E, as transcribed in the audit library moc4lib.CONTACT_SETS (same lexicographic 0-based permutation order, checked).",
    meaning="Each T_nu is a set of permutation indices at which the gap vector of Zheng's reduction vanishes (the permutations active on the supporting line at a contact). They are sets of permutations, not (a, b, U) configurations.",
    sets={str(k): sorted(v) for k, v in sorted(m4.CONTACT_SETS.items())}),
  numerical_contacts=num,
  numerical_contacts_note="Found by the constrained calibration run (dperm = dblock = 0.15 penalties); kept if dist_genperm > 0.1, dist_blockdiag > 0.1 and both spectra have separation > 0.1. All are VERTEX contacts: d(U) equals a hull vertex z_sigma at a U far from generalized permutations and block forms. The fiber of d over a vertex is positive-dimensional, so these are expected and are not evidence about MOC. No nontrivial EDGE contact was found.")
json.dump(out, open('contact_n4.json', 'w'), indent=1)
print(len(num), [(c['kind'], c['hp_m_60digits'][:8], c['hp_m_60digits'][-6:]) for c in num])
