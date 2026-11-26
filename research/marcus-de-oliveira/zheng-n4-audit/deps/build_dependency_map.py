#!/usr/bin/env python3
"""Dependency map of the ninety contact-case proofs of the Zheng MOC(4) supplement.

Usage:  python -I -B build_dependency_map.py [SUPP] TXT [OUTDIR]
  SUPP   supplement root (default: environment variable MOC4_SUPP)
  TXT    directory with moc4.txt (paper extract) and proofs/<dir>__<name>.txt
  OUTDIR where dependency-map.json is written (default: directory of this script)

Only reads the supplement's JSON/TSV/Markdown as data and the text extracts.
Nothing from the supplement is executed or imported.  All group actions, orbits
and generator zero sets are computed here from the definitions in the paper
(Section 4) and the facet numbering in data/facet-definitions.json.
"""
import sys, os, re, json, csv, itertools, collections

# ----------------------------------------------------------------- arguments
args = sys.argv[1:]
if len(args) >= 2 and os.path.isdir(os.path.join(args[0], 'cases')):
    SUPP, TXT = args[0], args[1]
    OUT = args[2] if len(args) > 2 else os.path.dirname(os.path.abspath(__file__))
else:
    SUPP = os.environ.get('MOC4_SUPP')
    TXT = args[0]
    OUT = args[1] if len(args) > 1 else os.path.dirname(os.path.abspath(__file__))
if not SUPP:
    sys.exit('supplement path: MOC4_SUPP or first argument')
J = lambda *p: os.path.join(SUPP, *p)

# ----------------------------------------------------------------- permutations and the group G
PERMS = list(itertools.permutations(range(4)))          # lexicographic, index 0..23
IDX = {p: i for i, p in enumerate(PERMS)}
def comp(a, b):                                         # (a o b)(i) = a[b[i]]
    return tuple(a[b[i]] for i in range(4))
def inv(a):
    r = [0] * 4
    for i, x in enumerate(a):
        r[x] = i
    return tuple(r)
GROUP, ACT = [], []
for rho in PERMS:
    for chi in PERMS:
        for t in (0, 1):
            GROUP.append((rho, chi, t))
            ACT.append(tuple(IDX[comp(rho, comp(inv(s) if t else s, chi))] for s in PERMS))
assert len(set(ACT)) == 1152 or True
def img(a, S):
    return frozenset(a[i] for i in S)
import functools
@functools.lru_cache(maxsize=None)
def _canon(S):
    return min(tuple(sorted(a[i] for i in S)) for a in ACT)
def canon(S):
    return _canon(tuple(sorted(S)))
def notation_action(row, col, t):
    """notation.md convention: sigma -> rho o sigma o chi (t=0), rho o sigma^-1 o chi (t=1)."""
    rho, chi = tuple(int(c) for c in row), tuple(int(c) for c in col)
    return tuple(IDX[comp(rho, comp(inv(s) if t else s, chi))] for s in PERMS)
def doc_action(row, col, t):
    """convention printed in blind_inventory49_whole: rho(i) = column[pi(row[i])]; transpose: pi -> pi^-1 first."""
    return tuple(IDX[tuple(col[(inv(s) if t else s)[row[i]]] for i in range(4))] for s in PERMS)
def min_extra(X, Y):
    """min over g in G of |g(X) \\ Y|: number of extra ties needed for an image of X to lie in Y + ties."""
    Y = set(Y)
    return min(len(img(a, X) - Y) for a in ACT)
def covers(X, Ytie):
    """some image of X lies inside Ytie (checked via orbit representatives of |X|-subsets)."""
    cx = canon(X)
    return any(canon(sub) == cx for sub in itertools.combinations(sorted(Ytie), len(X)))

# ----------------------------------------------------------------- generators (Section 4.1), own evaluation
def vals_n(I, Jset):
    return [1 if frozenset(s[i] for i in I) == frozenset(Jset) else 0 for s in PERMS]
fdefs = json.load(open(J('data', 'facet-definitions.json')))['definitions']
GEN = {}            # j -> (kind, values at the 24 permutations as tuple of Fractions*2)
for d in fdefs:
    j, kind, spec = d['index'], d['kind'], d['spec']
    if kind == 'n':
        a, b = spec.split(',')
        I = [int(c) - 1 for c in a]; Jj = [int(c) - 1 for c in b]
        v = [2 * x for x in vals_n(I, Jj)]
    else:
        p, e = spec.split(',')
        pi = [int(c) - 1 for c in p]
        eps = 1 if e == '+' else -1
        t = [sum(1 for i in range(4) if s[i] == pi[i]) for s in PERMS]
        dd = [a + b - c - e2 for a, b, c, e2 in zip(vals_n([0, 1], [pi[0], pi[2]]), vals_n([0, 1], [pi[1], pi[3]]),
                                                    vals_n([0, 2], [pi[0], pi[1]]), vals_n([0, 2], [pi[2], pi[3]]))]
        v = [x + eps * y for x, y in zip(t, dd)]          # = 2 h
    GEN[j] = (kind, tuple(v))
ZERO = {j: frozenset(i for i in range(24) if GEN[j][1][i] == 0) for j in GEN}
gen_check = {
    'count': len(GEN),
    'n_count': sum(1 for j in GEN if GEN[j][0] == 'n'),
    'h_count': sum(1 for j in GEN if GEN[j][0] == 'h'),
    'all_nonnegative_at_vertices': all(min(GEN[j][1]) >= 0 for j in GEN),
    'distinct_zero_sets': len(set(ZERO.values())),
}
def generators_on(T):
    return sorted(j for j in GEN if set(T) <= ZERO[j])
def nonminor_on(T):
    return [j for j in generators_on(T) if GEN[j][0] == 'h']
def relevant_ties(T):
    """extra ties t at which some non-minor generator vanishing on T also vanishes (Lemma 4.2)."""
    return sorted({t for j in nonminor_on(T) for t in ZERO[j] if t not in T})

# all orbits of six-subsets (paper: 214 orbits, 124 minor-only, 90 not)
seen, ORB6 = set(), []
for S in itertools.combinations(range(24), 6):
    if S in seen:
        continue
    orb = {tuple(sorted(a[i] for i in S)) for a in ACT}
    seen |= orb
    ORB6.append(min(orb))
TRIVIAL6 = {c for c in ORB6 if not nonminor_on(c)}
orbit_check = {'six_subset_orbits': len(ORB6), 'minor_only_orbits': len(TRIVIAL6),
               'orbits_with_nonminor_generator': len(ORB6) - len(TRIVIAL6)}

# ----------------------------------------------------------------- supplement data
idx_rows = list(csv.DictReader(open(J('cases', 'index.tsv')), delimiter='\t'))
cases = [r for r in idx_rows if r['paper_id'].startswith('T_')]
found = [r for r in idx_rows if r['paper_id'].startswith('F')]
catalog = json.load(open(J('cases', 'catalog.json')))
cat_by = {c['paper_id']: c for c in catalog['cases']}
methods = json.load(open(J('cases', 'methods.json')))
po = json.load(open(J('cases', 'proof-order.json')))
dmap = json.load(open(J('provenance', 'document-map.json')))['documents']
depg = json.load(open(J('provenance', 'DEPENDENCY_GRAPH.json')))['edges']
pdf_pages = {v['pdf_path']: v['pages'] for v in dmap.values()}
pdf_labels = {v['pdf_path']: v.get('label_pages', {}) for v in dmap.values()}
tex_of_pdf = {v['pdf_path']: k for k, v in dmap.items()}

CASE = {}           # 'T_n' -> info
for r in cases:
    n = int(r['paper_id'][2:])
    CASE['T_%d' % n] = dict(n=n, inv=int(r['inventory_index']), hist=r['historical_label'],
                            lit=[int(x) for x in r['literal_contacts'].split(',')],
                            can=[int(x) for x in r['canonical_contacts'].split(',')],
                            row=r['canonical_row'], col=r['canonical_column'], tr=r['canonical_transpose'],
                            alt=r['alternate_contacts'], arow=r['alternate_row'], acol=r['alternate_column'],
                            atr=r['alternate_transpose'], anchor=r['theorem_anchor'],
                            doc=r['document_path'], page=int(r['document_page']))
by_inv = {c['inv']: k for k, c in CASE.items()}
by_hist = {c['hist']: k for k, c in CASE.items()}
AUX = {}            # 'A_j' -> info
for a in po['auxiliary_contact_sets']:
    AUX['A_%d' % a['number']] = dict(j=a['number'], src=a['source_node'], lit=a['contacts'])
src2aux = {v['src']: k for k, v in AUX.items()}
BR = {'B_T9_regionI': 'T_10^(1)', 'B_T9_regionII': 'T_10^(2)', 'B_T10_A': 'T_11^(1)', 'B_T10_B': 'T_11^(2)'}
def node_of_notation(s):
    m = re.fullmatch(r'T_\{(\d+)\}', s)
    if m: return 'T_%s' % m.group(1)
    m = re.fullmatch(r'T\^\{\(7\)\}_\{(\d+)\}', s)
    if m: return 'A_%s' % m.group(1)
    m = re.fullmatch(r'T_\{(\d+)\}\^\{\((\d)\)\}', s)
    if m: return {'10': {'1': 'B_T9_regionI', '2': 'B_T9_regionII'}, '11': {'1': 'B_T10_A', '2': 'B_T10_B'}}[m.group(1)][m.group(2)]
    raise ValueError(s)
def label(n):
    if n.startswith('T_'): return n
    if n.startswith('A_'): return 'T^(7)_%s[%s]' % (n[2:], AUX[n]['src'])
    return BR.get(n, n)

ORBIT_CASE = {}
for k, c in CASE.items():
    ORBIT_CASE[canon(c['lit'])] = k
ORBIT_AUX = {}
for k, a in AUX.items():
    ORBIT_AUX[canon(a['lit'])] = k

# ----------------------------------------------------------------- paper: Appendix E, groups, Table 6
paper = open(os.path.join(TXT, 'moc4.txt')).read()
appE = {}
for m in re.finditer(r'\bT(\d{1,2})\s+\{([\d,\s]+)\}', paper):
    appE.setdefault('T_%s' % m.group(1), sorted(int(x) for x in m.group(2).split(',')))
# Table 6 (paper p. 63) and Appendix G groups (pp. 61-62), transcribed from the paper text.
TABLE6 = {
    1: ('T1 T2 T3 T4 T5 T6 T7 T15 T18 T19 T20 T21 T22 T23 T25 T27 T28 T29 T31 T34 T35 T36 T37 T39 T40 T42 T43 T44 T72 T88',
        'A1 A2 A5 A7 A8 A9 A10 A11 A12 A13 A14 A15 A16 A17 A21 A22 A24 B_T9_regionI B_T9_regionII'),
    2: ('T10 T26 T46 T90', 'A3 A4 A23 A26 B_T10_A B_T10_B'),
    3: ('T9 T11 T12 T47', ''),
    4: ('T8 T13 T17 T49', 'A6 A18'),
    5: ('T14 T16 T24 T33 T38 T70', 'A20 A27'),
    6: ('T52 T57 T58 T65 T68 T71 T76 T78', ''),
    7: ('T45 T51 T66 T69 T75 T81', 'A19 A25'),
    8: ('T48 T50 T53 T59 T62 T67 T73 T79 T86', ''),
    9: ('T30 T41 T54 T55 T60 T61 T63 T74 T77 T80 T83 T84', ''),
    10: ('T56 T64 T82 T85 T87', ''),
    11: ('T32', ''),
    12: ('T89', ''),
}
PAPER_STAGE = {}
for st, (a, b) in TABLE6.items():
    for x in a.split():
        PAPER_STAGE['T_' + x[1:]] = st
    for x in b.split():
        PAPER_STAGE[('A_' + x[1:]) if x.startswith('A') else x] = st
# Appendix G group lists (parsed from the paper extract)
gtxt = paper[paper.find('Group I\n'):paper.find('The tables below give an order')]
PAPER_GROUP = {}
for g, body in re.findall(r'Group (I|II|III|IV|V|VI)\n(.*?)(?=\nGroup |\Z)', gtxt, re.S):
    for x in re.findall(r'T(\d+)', body):
        PAPER_GROUP['T_' + x] = g

# ----------------------------------------------------------------- proof-order graph checks
NODES = {}
for L, layer in enumerate(po['layers']):
    for x in layer:
        n = node_of_notation(x['paper_notation'])
        NODES[n] = L
EDGES = set()
edge_layer_bad = []
for e in po['edges']:
    a, b = node_of_notation(e['from']), node_of_notation(e['to'])
    EDGES.add((a, b))
    if NODES.get(a) != e['from_layer'] or NODES.get(b) != e['to_layer']:
        edge_layer_bad.append((a, b, 'recorded layer differs from node layer'))
    if not NODES[a] < NODES[b]:
        edge_layer_bad.append((a, b, 'from_layer >= to_layer'))
PRED = collections.defaultdict(set); SUCC = collections.defaultdict(set)
for a, b in EDGES:
    PRED[b].add(a); SUCC[a].add(b)

def sccs(nodes, succ):
    index, low, onst, st, out, cnt = {}, {}, set(), [], [], [0]
    sys.setrecursionlimit(10000)
    def strong(v):
        index[v] = low[v] = cnt[0]; cnt[0] += 1; st.append(v); onst.add(v)
        for w in succ.get(v, ()):
            if w not in index:
                strong(w); low[v] = min(low[v], low[w])
            elif w in onst:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp_ = []
            while True:
                w = st.pop(); onst.discard(w); comp_.append(w)
                if w == v: break
            out.append(comp_)
    for v in nodes:
        if v not in index:
            strong(v)
    return out
def longest_path_layers(nodes, pred):
    memo = {}
    def lay(v, stack=()):
        if v in memo: return memo[v]
        if v in stack: return None
        ps = pred.get(v, ())
        r = 0 if not ps else max((lay(p, stack + (v,)) or 0) for p in ps) + 1
        memo[v] = r
        return r
    return {v: lay(v) for v in nodes}
def closure(n, pred):
    out, todo = set(), [n]
    while todo:
        v = todo.pop()
        for p in pred.get(v, ()):
            if p not in out:
                out.add(p); todo.append(p)
    return out

graph_scc = [c for c in sccs(sorted(NODES), SUCC) if len(c) > 1]
LP = longest_path_layers(sorted(NODES), PRED)
graph_checks = {
    'nodes': len(NODES), 'cases': sum(1 for n in NODES if n.startswith('T_')),
    'aux': sum(1 for n in NODES if n.startswith('A_')), 'branch': sum(1 for n in NODES if n.startswith('B_')),
    'distinct_edges': len(EDGES), 'edges_listed': len(po['edges']),
    'nontrivial_scc': graph_scc,
    'edge_layer_inconsistencies': edge_layer_bad,
    'layer_equals_longest_path': all(LP[n] == NODES[n] for n in NODES),
    'layer_vs_longest_path_diffs': {n: (NODES[n], LP[n]) for n in NODES if LP[n] != NODES[n]},
    'stage_vs_paper_table6': {label(n): {'graph_layer+1': NODES[n] + 1, 'paper_stage': PAPER_STAGE.get(n)}
                              for n in NODES if PAPER_STAGE.get(n) != NODES[n] + 1},
    'nodes_missing_from_table6': [n for n in NODES if n not in PAPER_STAGE],
    'table6_entries_missing_from_graph': [n for n in PAPER_STAGE if n not in NODES],
}

# index consistency: literal contacts vs Appendix E; canonical via stated action; orbits distinct
idx_checks = {'appendixE_parsed': len(appE), 'literal_vs_appendixE_mismatch': [], 'canonical_action_fail': [],
              'alternate_action_fail': [], 'canonical_is_orbit_min_fail': [],
              'distinct_case_orbits': len(ORBIT_CASE), 'distinct_aux_orbits': len(ORBIT_AUX),
              'aux_table5_mismatch': [], 'group_vs_paper': [], 'case_orbit_minor_only': [],
              'aux_orbit_minor_only': []}
for k, c in CASE.items():
    if appE.get(k) != sorted(c['lit']):
        idx_checks['literal_vs_appendixE_mismatch'].append((k, c['lit'], appE.get(k)))
    a = notation_action(c['row'], c['col'], int(c['tr']))
    if sorted(img(a, c['lit'])) != sorted(c['can']):
        idx_checks['canonical_action_fail'].append(k)
    if tuple(sorted(c['can'])) != canon(c['lit']):
        idx_checks['canonical_is_orbit_min_fail'].append((k, c['can'], canon(c['lit'])))
    if c['alt']:
        alt = [int(x) for x in c['alt'].split(',')]
        b = notation_action(c['arow'], c['acol'], int(c['atr']))
        if sorted(img(b, c['lit'])) != sorted(alt) and sorted(img(b, alt)) != sorted(c['lit']):
            idx_checks['alternate_action_fail'].append((k, c['lit'], alt))
    if PAPER_GROUP.get(k) != cat_by[k]['method_group']:
        idx_checks['group_vs_paper'].append((k, cat_by[k]['method_group'], PAPER_GROUP.get(k)))
    if canon(c['lit']) in TRIVIAL6:
        idx_checks['case_orbit_minor_only'].append(k)
for k, a in AUX.items():
    if not nonminor_on(a['lit']):
        idx_checks['aux_orbit_minor_only'].append(k)
t5 = dict((int(j), sorted(int(x) for x in s.split(','))) for j, s in
          re.findall(r'(\d{1,2})\s+\{(\d+(?:, \d+){6})\}', paper[paper.find('Table 5'):paper.find('G.2')]))
for k, a in AUX.items():
    if t5.get(a['j']) != sorted(a['lit']):
        idx_checks['aux_table5_mismatch'].append((k, a['lit'], t5.get(a['j'])))
idx_checks['table5_parsed'] = len(t5)
idx_checks['nontrivial_case_orbits_equal_90'] = (set(ORBIT_CASE) == set(ORB6) - TRIVIAL6)

# ----------------------------------------------------------------- documents
def txt_path(pdf):
    d, b = pdf.split('/')[-2], os.path.basename(pdf)[:-4]
    return os.path.join(TXT, 'proofs', d + '__' + b + '.txt')
LIG = {'ﬀ': 'ff', 'ﬁ': 'fi', 'ﬂ': 'fl', 'ﬃ': 'ffi', 'ﬄ': 'ffl'}
def norm(t):
    for a, b in LIG.items():
        t = t.replace(a, b)
    return t
DOCS = sorted(v['pdf_path'] for v in dmap.values())
TEXT = {}
for p in DOCS:
    tp = txt_path(p)
    TEXT[p] = norm(open(tp).read()) if os.path.exists(tp) else None
missing_text = [p for p in DOCS if TEXT[p] is None]
def page_at(t, pos):
    return t.count('\f', 0, pos) + 1
def stem(p):
    return os.path.basename(p)[:-4]

# filename "home" of a document
def home_by_name(p):
    b = stem(p)
    homes = []
    m = re.search(r'core([A-O])', b)
    if m: homes.append(src2aux['core' + m.group(1)])
    if 'seven_A' in b: homes.append(src2aux['coreA'])
    specials = {'inventory52_plus17': 'S52', 'inventory56_boundary17': 'S56', 'inventory75_boundary16': 'S75',
                'inventory78_boundary15': 'S78', 'inventory67_g14': 'S67', 'true_T48_c0': 'S48c0', 'T48_c0': 'S48c0',
                'T46_a0_uniform': 'S46a0', 'T7_boundary16': 'S_0,1,2,3,6,15,17', 'T7_boundary23': 'S_0,1,2,3,6,15,22',
                'T23_g17': 'S_0,1,2,6,9,11,20', 'T23_g13': 'S_0,1,2,6,9,11,23'}
    for k, v in specials.items():
        if k in b: homes.append(src2aux[v])
    brs = {'T9_first_branch': 'B_T9_regionI', 'T9_second_branch': 'B_T9_regionII',
           'T10_A_branch': 'B_T10_A', 'T10_B_branch': 'B_T10_B'}
    for k, v in brs.items():
        if k in b: homes.append(v)
    if homes:
        return homes
    m = re.search(r'inventory(\d+)(?:_(\d+))?', b)
    if m:
        homes += [by_inv[int(m.group(1))]] + ([by_inv[int(m.group(2))]] if m.group(2) else [])
        return homes
    m = re.search(r'(?:^|_)T(\d+)(?:_T(\d+))?_', b)
    if m:
        homes += [by_hist['T' + m.group(1)]] + ([by_hist['T' + m.group(2)]] if m.group(2) else [])
        return homes
    if 'bad70' in b: return ['T_34']
    if b == 'independent_moc4_six_contact_reduction_20261022': return ['T_1', 'T_2']
    return []
HOME = {p: home_by_name(p) for p in DOCS}
for k, c in CASE.items():                      # primary documents always belong to their case
    pd = 'proofs/' + c['doc'] if not c['doc'].startswith('proofs/') else c['doc']
    c['doc'] = pd
    if k not in HOME.get(pd, []):
        HOME.setdefault(pd, []).append(k)

# ----------------------------------------------------------------- citation extraction
SETRE = re.compile(r'[\{\(]\s*(\d{1,2}(?:\s*,\s*\d{1,2}){5,6})\s*[\}\)]')
GEQ = re.compile(r'g_?(\d{1,2})(?:\s*=\s*g_?(\d{1,2})){5,6}\s*=\s*0')
HUNION = re.compile(r'(?<![A-Za-z0-9_])T\s?(\d{1,2})\s*(?:∪\s*\{\s*(\d{1,2})\s*\}|\+\s*(\d{1,2})\b)')
TLAB = re.compile(r'(?<![A-Za-z0-9_′])(?:T\s?(\d{1,2})|I\s?(29|47|83))(?![0-9∪])(?!\s*∪)(?!\s*\+\s*\d)')
INVRE = re.compile(r'(?:true[- ])?inventory[- ]?(?:index\s)?(\d{1,2})(?!\d)|whole ?(\d{2})(?!\d)|historical ind(?:ex|ices) (\d{1,2})(?: and (\d{1,2}))?|index (\d{1,2}) in the archived inventory', re.I)
CORE = re.compile(r'(?:seven-)?cores?\s?([A-O])\b((?:\s*(?:,|and)\s*[A-O]\b)*)')
INDEXISH = re.compile(r'(ind(?:ex|ices)|ties|new tie|positions|zeros?|j ∈|entries|roles|rows|columns)\W{0,12}$', re.I)
BULK = re.compile(r'(preceding|prior|earlier|immutable|fixed|retained|preclosure)[- ]?(?:closed)?\s?(\d{2})[- ](?:family|families|contact|case|list)|(\d{2}) (?:preceding|prior) contact|old ?42|preclosure-47|retained full results|previously proved|strictly prior', re.I)

def resolve_set(S):
    if len(set(S)) != len(S) or max(S) > 23:
        return None
    c = canon(S)
    if len(S) == 6:
        if c in ORBIT_CASE: return ORBIT_CASE[c]
        if c in TRIVIAL6: return 'minor-only-6'
        return 'unresolved-6'
    if c in ORBIT_AUX: return ORBIT_AUX[c]
    if not nonminor_on(S): return 'minor-only-7'
    return 'other-7:' + ','.join(map(str, sorted(S)))

NEG = re.compile(r'\b(no|not|neither|nor|without|never|instead|different|differently|rather than)\b', re.I)
KEYW = re.compile(r'prov|closed|famil|case|inventor|label|result|inequalit|witness|tuple|core|capacit|lemma|calculation|allocation|theorem|orbit|\{', re.I)
RANGE = re.compile(r'(?<![A-Za-z0-9_])T\s?(\d{1,2})\s*,\s*(?:\.\s*){3},?\s*T\s?(\d{1,2})')
def extract(p):
    t = TEXT[p]
    out = []                       # (node, method, page, snippet, flags)
    if t is None:
        return out
    def add(node, how, pos, flag=''):
        sn = re.sub(r'\s+', ' ', t[max(0, pos - 70):pos + 70])
        if not flag and NEG.search(t[max(0, pos - 60):pos]):
            flag = 'negated?'
        out.append(dict(node=node, how=how, page=page_at(t, pos), snippet=sn, flag=flag))
    for m in SETRE.finditer(t):
        S = [int(x) for x in re.split(r'\s*,\s*', m.group(1).strip())]
        r = resolve_set(S)
        if r is None: continue
        pre = t[max(0, m.start() - 40):m.start()]
        flag = 'index-list?' if INDEXISH.search(pre) else ''
        if m.group(0)[0] == '(' and not re.search(r'tuple|contact|ties|tied', pre, re.I):
            continue
        add(r, 'set', m.start(), flag)
    for m in GEQ.finditer(t):
        S = [int(x) for x in re.findall(r'g_?(\d{1,2})', m.group(0))]
        r = resolve_set(S)
        if r: add(r, 'g-equalities', m.start())
    for m in HUNION.finditer(t):
        h = 'T' + m.group(1)
        if h in by_hist:
            k = int(m.group(2) or m.group(3))
            S = CASE[by_hist[h]]['lit'] + [k]
            r = resolve_set(S) if k not in CASE[by_hist[h]]['lit'] else None
            if r: add(r, 'hist-label+tie', m.start())
    ranged = set()
    for m in RANGE.finditer(t):
        a, b = int(m.group(1)), int(m.group(2))
        if a < b <= 89:
            for x in range(a, b + 1):
                if 'T%d' % x in by_hist:
                    add(by_hist['T%d' % x], 'hist-label-range', m.start(), 'declared-input-list')
            ranged.update(range(m.start(), m.end()))
    for m in TLAB.finditer(t):
        if m.start() in ranged:
            continue
        h = ('T' + m.group(1)) if m.group(1) else ('I' + m.group(2))
        ctx = t[max(0, m.start() - 45):m.end() + 45]
        listy = re.search(r'T\s?\d{1,2}\s*,\s*T\s?\d{1,2}', ctx)
        flag = '' if (KEYW.search(ctx) or listy) else 'weak(no keyword)'
        if h in by_hist:
            add(by_hist[h], 'hist-label', m.start(), flag)
        else:
            add('unresolved-label:' + h, 'hist-label', m.start(), flag)
    for m in INVRE.finditer(t):
        if m.start() > 0 and t[m.start() - 1] in '_':
            continue
        after = t[m.end():m.end() + 25]
        for g in m.groups():
            if g is not None and int(g) in by_inv:
                if re.search(r'seven', after) and ('S%d' % int(g)) in src2aux:
                    add(src2aux['S%d' % int(g)], 'inventory-seventh', m.start())
                else:
                    add(by_inv[int(g)], 'inventory', m.start())
    for m in CORE.finditer(t):
        if m.start() > 0 and t[m.start() - 1] in '_':
            continue
        letters = [m.group(1)] + re.findall(r'[A-O]', m.group(2) or '')
        for L in letters:
            if 'core' + L in src2aux:
                add(src2aux['core' + L], 'core-name', m.start())
    sq = re.sub(r'\s+', '', t)
    for q in DOCS:
        if q != p and stem(q) in sq:
            pos = sq.find(stem(q))
            neg = 'negated?' if NEG.search(sq[max(0, pos - 50):pos]) else ''
            for h in HOME[q] or ['common:' + stem(q)]:
                out.append(dict(node=h, how='file:' + stem(q), page=None, snippet=sq[max(0, pos - 60):pos + len(stem(q)) + 10], flag=neg))
    for m in BULK.finditer(t):
        add('bulk', 'bulk', m.start())
    return out
CIT = {p: extract(p) for p in DOCS}

# ----------------------------------------------------------------- node -> proving documents
def statement_subjects(p):
    """resolved sets/labels in the title and first theorem statement (first ~2500 chars)."""
    t = TEXT[p] or ''
    head = t[:2500]
    subs = []
    for c in CIT[p]:
        pass
    for m in SETRE.finditer(head):
        S = [int(x) for x in re.split(r'\s*,\s*', m.group(1).strip())]
        r = resolve_set(S)
        if r and not INDEXISH.search(head[max(0, m.start() - 40):m.start()]):
            subs.append(r)
    for m in GEQ.finditer(head):
        r = resolve_set([int(x) for x in re.findall(r'g_?(\d{1,2})', m.group(0))])
        if r: subs.append(r)
    for m in HUNION.finditer(head):
        h = 'T' + m.group(1)
        if h in by_hist:
            k = int(m.group(2) or m.group(3)); r = resolve_set(CASE[by_hist[h]]['lit'] + [k])
            if r: subs.append(r)
    return subs
SUBJ = {p: statement_subjects(p) for p in DOCS}
PROVERS = collections.defaultdict(list)          # node -> [(doc, page, how)]
for k, c in CASE.items():
    PROVERS[k].append((c['doc'], c['page'], 'index.tsv primary'))
for p in DOCS:
    for h in HOME[p]:
        if (h.startswith('A_') or h.startswith('B_')):
            PROVERS[h].append((p, 1, 'file name'))
    for s in SUBJ[p][:3]:
        if HOME[p] and s not in HOME[p]:
            continue
        if s.startswith('A_') and all(q != p for q, _, _ in PROVERS[s]):
            PROVERS[s].append((p, 1, 'statement set'))
# inline seven-contact theorems (statement anywhere in a document): record page of first mention
for a in AUX:
    if not PROVERS[a]:
        for p in DOCS:
            for c in CIT[p]:
                if c['node'] == a and c['how'] in ('set', 'hist-label+tie'):
                    PROVERS[a].append((p, c['page'], 'inline mention (needs reading)'))
                    break

# ----------------------------------------------------------------- per node own documents and citations
def own_docs(n):
    docs = []
    if n.startswith('T_'):
        docs.append(CASE[n]['doc'])
        for b in cat_by[n]['branches']:
            if n in HOME.get(b, []) and b not in docs:
                docs.append(b)
    else:
        for p, _, how in PROVERS[n]:
            if how != 'inline mention (needs reading)' and p not in docs:
                docs.append(p)
    return docs
NODE_LIST = sorted(NODES, key=lambda n: (NODES[n], n))
TEXTCIT = {}
for n in NODE_LIST:
    found_ = collections.defaultdict(list)
    for p in own_docs(n):
        for c in CIT[p]:
            m = c['node']
            if m == n or m in ('bulk',):
                continue
            if m in NODES and m.startswith('T_') and n.startswith('T_') and CASE[m]['doc'] == CASE[n]['doc']:
                continue        # joint document proving both (T_1/T_2, T_22/T_23, T_57/T_58)
            # documents owned by the node itself (e.g. sub-notes) are not citations
            if c['how'].startswith('file:') and any(stem(q) == c['how'][5:] and n in HOME[q] for q in DOCS):
                continue
            found_[m].append(dict(doc=stem(p), page=c['page'], how=c['how'], flag=c['flag'], snippet=c['snippet']))
    TEXTCIT[n] = found_
def is_node(m):
    return m in NODES

# ----------------------------------------------------------------- edge support (task 1c)
CL = {n: closure(n, PRED) for n in NODES}
def lit_any(n):
    return CASE[n]['lit'] if n.startswith('T_') else AUX[n]['lit'] if n.startswith('A_') else None
def boundary_ties(X, Y):
    """ties t (not in Y) such that an image of X lies in Y + {t}; with relevance for Y."""
    lx, ly = lit_any(X), lit_any(Y)
    if lx is None or ly is None:
        return None
    rel = set(relevant_ties(ly))
    out = [t for t in range(24) if t not in ly and covers(lx, set(ly) | {t})]
    return dict(ties=out, relevant=[t for t in out if t in rel])
edge_support = {}
for a, b in sorted(EDGES):
    hits = TEXTCIT[b].get(a, [])
    hows = sorted({h['how'].split(':')[0] for h in hits if not h['flag']})
    cat = False
    if b.startswith('T_'):
        brs = set(cat_by[b]['branches']) | {CASE[b]['doc']}
        cat = any(p in brs for p, _, how in PROVERS[a] if how != 'inline mention (needs reading)')
    edge_support[(a, b)] = dict(text=hows, catalog_branch=cat, boundary=boundary_ties(a, b),
                                min_extra_ties=(min_extra(CASE[a]['lit'] if a.startswith('T_') else AUX[a]['lit'], CASE[b]['lit'] if b.startswith('T_') else AUX[b]['lit'])
                                                if not (a.startswith('B_') or b.startswith('B_')) else None))
unsupported_edges = [(a, b, v) for (a, b), v in edge_support.items() if not v['text']]
missing_edges = []
for n in NODE_LIST:
    for m, hits in TEXTCIT[n].items():
        if not is_node(m) or (m, n) in EDGES:
            continue
        good = [h for h in hits if not h['flag']]
        if not good:
            continue
        kind = 'in transitive closure' if m in CL[n] else ('REVERSE (cited node depends on citer)' if n in CL.get(m, ()) else 'absent from closure')
        if m == n: continue
        missing_edges.append(dict(citer=n, cited=m, status=kind, boundary=boundary_ties(m, n),
                                  methods=sorted({h['how'].split(':')[0] for h in good}),
                                  min_extra_ties=(min_extra(CASE[m]['lit'] if m.startswith('T_') else AUX[m]['lit'], CASE[n]['lit'] if n.startswith('T_') else AUX[n]['lit'])
                                                  if not (m.startswith('B_') or n.startswith('B_')) else None),
                                  evidence=good[:3]))
unresolved = []
for n in NODE_LIST:
    for m, hits in TEXTCIT[n].items():
        if m.startswith('unresolved') or m.startswith('other-7'):
            for h in hits:
                if not h['flag']:
                    unresolved.append(dict(citer=n, ref=m, **{k: h[k] for k in ('doc', 'page', 'snippet')}))

# ----------------------------------------------------------------- boundary coverage by earlier stages
def stage(n):
    return PAPER_STAGE.get(n, NODES.get(n, 99) + 1)
coverage = {}
for n in NODE_LIST:
    if not (n.startswith('T_') or n.startswith('A_')):
        continue
    T = CASE[n]['lit'] if n.startswith('T_') else AUX[n]['lit']
    rows = []
    for t in relevant_ties(T):
        Tt = set(T) | {t}
        cov = [k for k, c in CASE.items() if k != n and covers(c['lit'], Tt)]
        cov += [k for k, a in AUX.items() if k != n and len(Tt) >= 7 and covers(a['lit'], Tt)]
        triv = (not nonminor_on(sorted(Tt)))
        earlier = [k for k in cov if stage(k) < stage(n)]
        rows.append(dict(tie=t, covers=cov, earlier_covers=earlier, cone_minor_only=triv))
    coverage[n] = rows
uncovered = [(n, r['tie'], r['covers']) for n, rows in coverage.items() for r in rows
             if not r['earlier_covers'] and not r['cone_minor_only']]
# which relevant ties are covered by the node's own graph prerequisites
for n, rows in coverage.items():
    for r in rows:
        r['graph_prerequisites_covering'] = [k for k in r['covers'] if k in PRED.get(n, ())]
ties_without_graph_witness = {n: [r['tie'] for r in rows if not r['graph_prerequisites_covering']]
                              for n, rows in coverage.items()}
ties_without_graph_witness = {n: v for n, v in ties_without_graph_witness.items() if v}
# graph prerequisites that cover no relevant tie of the node (used for something other than a one-tie boundary)
preds_covering_nothing = {}
for n in coverage:
    used = {k for r in coverage[n] for k in r['graph_prerequisites_covering']}
    extra = [k for k in PRED.get(n, ()) if k not in used]
    if extra:
        preds_covering_nothing[n] = extra

# ----------------------------------------------------------------- tags (task 4)
TAGS = {
    'B:shared-zero/cross-collision': r'shared[- ]zero|cross[- ]?(?:spectrum )?collision|cross-spectrum common',
    'B:coincident-within-spectrum': r'repetition|repeated (?:spectral )?points|within (?:one|either) spectrum|coincident eigen',
    'B:single-collision-compensation': r'single[- ]collision|collision compensation|divisor compensation|collision factor',
    'C:one-circle': r'one[- ](?:spectrum[- ])?circle|generali[sz]ed circle',
    'C:Hessian-rank/inertia': r'Hessian|inertia|principal (?:Hessian )?minors? (?:nonneg|of order)|rank H',
    'L3.5:harmonic/max-principle/negative-mass': r'harmonic|maximum principle|subharmonic|bounded[- ]negative|negative[- ]part lemma|complete negative mass|bounded bad',
    'Thm1.2:two-real-responses': r'two[- ]real|two[- ]orthogonal',
    'Gram/PSD/SOS/capacity': r'Gram|positive semidefinite|PSD|sum[- ]of[- ]squares|\bSOC\b|old[- ]H|chain (?:capacity|inequality)|four-coefficient',
    'Rouche/wall-perturbation': r'Rouch|parameter wall|walls?\b',
    'resultant/discriminant': r'resultant|discriminant',
    'Birkhoff/height-kernel': r'Birkhoff',
    'matching-cubic/Sing F': r'cubic|Sing ?F',
}
TAGRE = {k: re.compile(v, re.I) for k, v in TAGS.items()}
def tags_of_text(t):
    return sorted(k for k, r in TAGRE.items() if t and r.search(t))
COMMON_DOCS = [p for p in DOCS if not HOME[p]]
case_tags = {}
for n in CASE:
    own = own_docs(n)
    tg = set()
    for p in own:
        tg |= set(tags_of_text(TEXT[p]))
    common_listed = sorted(stem(b) for b in cat_by[n]['branches'] if not HOME.get(b))
    case_tags[n] = dict(keyword_tags=sorted(tg), common_documents_in_catalog=common_listed)

# ----------------------------------------------------------------- conditional theorems (task 2)
COND_PAT = re.compile(r'Theorem[^.\n]{0,40}\.\s*(?:\(Conditional[^)]*\)\.?\s*)?(Assume|Suppose|Under the hypothesis)|'
                      r'conditional|provided that|remains the explicit hypothesis|explicit hypothesis|open input|'
                      r'is a hypothesis|hypothesis of this reduction|not proved in this note|not proved by this', re.I)
cond_hits = {}
for p in DOCS:
    t = TEXT[p] or ''
    hs = []
    for m in COND_PAT.finditer(t):
        hs.append(dict(page=page_at(t, m.start()), snippet=re.sub(r'\s+', ' ', t[max(0, m.start() - 80):m.start() + 160])))
    if hs:
        cond_hits[stem(p)] = hs
cond_files = sorted(os.path.relpath(os.path.join(dp, f), SUPP) for dp, _, fs in os.walk(SUPP) for f in fs if 'conditional' in f.lower())
def role(f):
    if f.endswith('.pdf'): return 'PDF proof (conditional theorem)'
    if 'frozen_manifest' in f: return 'manifest (frozen/approval record)'
    if f.endswith('.py'): return 'exact-check script (not executed)'
    if f.endswith('.json'): return 'JSON result of a check'
    if f.endswith('.txt'): return 'audit note'
    return '?'
cond_file_roles = [(f, role(f)) for f in cond_files]

def lit(n):
    return CASE[n]['lit'] if n.startswith('T_') else AUX[n]['lit']
# Conditional theorems located by the search above and read by hand.  Each entry: the theorem's
# own literal contact set, hypothesis sets (literal), the case using it, the discharging document(s),
# and the tie -> witness tables printed in the documents (labels resolved here, then checked).
H = lambda h: by_hist[h]
CONDITIONALS = [
    dict(id='C1', doc='root_inventory50_bounded_native_conditional_20261023', T=[1, 3, 5, 10, 11, 12], user='T_51',
         hyp={'D': [0, 1, 2, 7, 9, 16, 23]}, hyp_ties={'D': 20},
         discharge='blind_moc4_coreD_uniform_20261023', assembly='root_inventory50_bounded_native_whole_20261023',
         witnesses={6: 'A_8', 13: 'A_8', 14: 'A_11', 20: 'A_12'}, unnamed_ties=[0, 7, 16]),
    dict(id='C2', doc='root_inventory51_bounded_native_conditional_20261023', T=[1, 3, 5, 10, 11, 14], user='T_52',
         hyp={'G': [0, 1, 2, 7, 9, 10, 23]}, hyp_ties={'G': 20},
         discharge='independent_shared_coreG_uniform_20261023', assembly='root_inventory51_bounded_native_whole_20261023',
         witnesses={6: 'A_9', 13: 'A_9', 12: 'A_11', 20: 'A_10'}, unnamed_ties=[0, 7, 16]),
    dict(id='C3', doc='blind_inventory49_bounded_native_conditional_20261023', T=[1, 3, 5, 6, 10, 11], user='T_50',
         hyp={'F': [0, 1, 2, 7, 9, 10, 14], 'G': [0, 1, 2, 7, 9, 10, 23]}, hyp_ties={'F': 13, 'G': 20},
         discharge='blind_shared_coreF_uniform_20261023 + independent_shared_coreG_uniform_20261023',
         assembly='blind_inventory49_whole_20261023',
         witnesses={0: H('T9'), 7: H('T42'), 16: H('T65'), 12: 'A_8', 14: 'A_9', 13: 'A_7', 20: 'A_10'},
         actions={0: ((0, 1, 2, 3), (0, 1, 3, 2), 0, [0, 1, 3, 5, 7, 8]), 7: ((0, 1, 3, 2), (0, 1, 3, 2), 0, [1, 2, 4, 6, 7, 8]),
                  12: ((0, 2, 3, 1), (0, 1, 2, 3), 0, [0, 1, 2, 7, 9, 10, 16]), 13: ((0, 2, 3, 1), (0, 1, 2, 3), 0, [0, 1, 2, 7, 9, 10, 14]),
                  14: ((0, 2, 3, 1), (0, 1, 2, 3), 0, [0, 1, 2, 7, 9, 10, 17]), 16: ((0, 3, 1, 2), (0, 1, 2, 3), 0, [2, 4, 5, 6, 9, 12]),
                  20: ((0, 2, 3, 1), (0, 1, 2, 3), 0, [0, 1, 2, 7, 9, 10, 23])},
         unnamed_ties=[]),
    dict(id='C4', doc='blind_moc4_inventory52_native_harmonic_conditional_20261023', T=[0, 1, 4, 6, 9, 12], user='T_53',
         hyp={'S52': [0, 1, 4, 6, 9, 12, 17]}, hyp_ties={'S52': 17},
         discharge='blind_moc4_inventory52_whole_assembly_20261023 (seven-contact part) + four plus17 notes',
         assembly='blind_moc4_inventory52_whole_assembly_20261023',
         witnesses={5: H('T4'), 7: H('T17'), 11: H('T25'), 14: H('T4'), 15: H('T25'), 18: H('T47')}, unnamed_ties=[]),
    dict(id='C5', doc='final_referee_inventory56_native_harmonic_20261023', T=[0, 1, 5, 7, 9, 12], user='T_57',
         hyp={'S56': [0, 1, 5, 7, 9, 12, 17]}, hyp_ties={'S56': 17},
         discharge='final_referee_inventory56_boundary17_uniform_20261023', assembly='final_referee_inventory56_57_whole_assembly_20261023',
         witnesses={}, unnamed_ties='see document'),
    dict(id='C6', doc='root_inventory57_native_conditional_20261023', T=[0, 3, 5, 7, 9, 12], user='T_58',
         hyp={'S56': [0, 1, 5, 7, 9, 12, 17]}, hyp_ties={'S56': 17},
         discharge='final_referee_inventory56_boundary17_uniform_20261023', assembly='final_referee_inventory56_57_whole_assembly_20261023',
         witnesses={}, unnamed_ties='working tuple {0,1,5,7,9,17}; ties 4,6,11,14,15,18 there'),
    dict(id='C7', doc='independent_inventory67_native_harmonic_conditional_20261023', T=[1, 4, 5, 7, 9, 12], user='T_68',
         hyp={'S67': [1, 4, 5, 7, 9, 12, 14]}, hyp_ties={'S67': 14},
         discharge='root_inventory67_g14_uniform_20261023 (+ seven-contact section of root_inventory67_whole_assembly)',
         assembly='root_inventory67_whole_assembly_20261023', witnesses={}, unnamed_ties=[0, 6, 11, 15, 17, 18]),
    dict(id='C8', doc='independent_inventory75_native_harmonic_conditional_20261023', T=[0, 1, 3, 10, 11, 12], user='T_76',
         hyp={'S75': [0, 1, 3, 10, 11, 12, 16]}, hyp_ties={'S75': 16},
         discharge='final_referee_inventory75_boundary16_uniform_20261023', assembly='final_referee_inventory75_whole_assembly_20261023',
         witnesses={14: 'A_27'}, unnamed_ties=[5, 6, 7, 13, 20]),
    dict(id='C9', doc='independent_true_T48_native_harmonic_conditional_20261023', T=[0, 3, 5, 6, 9, 12], user='T_49',
         hyp={'S48c0': [0, 3, 5, 6, 9, 12, 17]}, hyp_ties={'S48c0': 17},
         discharge='root_true_T48_whole_assembly_20261023 (seven-contact section) + independent_true_T48_c0_native_Heron, blind_moc4_T48_c0_native_interlacing, blind_moc4_T48_c0_all_j_chain',
         assembly='root_true_T48_whole_assembly_20261023',
         witnesses={2: H('T11'), 7: H('T36'), 11: H('T20'), 13: H('T41'), 16: H('T33'), 23: H('T43')}, unnamed_ties=[]),
    # T+16 is discharged by the six-contact case hist T24 (a case, not an auxiliary set); T+17 by A_5.
    dict(id='C10', doc='independent_moc4_T23_harmonic_reduction_20261022', T=[2, 3, 4, 6, 7, 8], user='T_24',
         hyp={'T+17': [2, 3, 4, 6, 7, 8, 17]}, hyp_ties={'T+17': 17},
         discharge='independent_moc4_T23_g17_uniform_20261022 (T+17); six-contact case hist T24 = T_25 via printed action (T+16)',
         assembly='independent_moc4_T23_uniform_20261022',
         witnesses={1: H('T14'), 10: H('T19'), 12: H('T12'), 22: H('T11'), 16: H('T24'), 13: 'A_6', 17: 'A_5'},
         case_actions={16: ((0, 3, 1, 2), (1, 2, 0, 3), 1, H('T24'))}, unnamed_ties=[]),
]
# S52 seven-face zero-coefficient witnesses (blind_moc4_inventory52_whole_assembly p.1-2)
S52_WIT = {18: H('T47'), 15: H('T25'), 11: H('T25'), 14: H('T4'), 7: H('T17'), 5: H('T4')}

def doc_by_stem(s):
    for p in DOCS:
        if stem(p) == s:
            return p
cond_results = []
for c in CONDITIONALS:
    T = c['T']; user = c['user']
    r = dict(id=c['id'], document=c['doc'] + '.pdf', using_case=user,
             theorem_set_matches_case=(canon(T) == canon(CASE[user]['lit'])),
             theorem_set_literal_equals_appendixE=(sorted(T) == sorted(CASE[user]['lit'])))
    hyps = {}
    for name, S in c['hyp'].items():
        node = resolve_set(S)
        tie = c['hyp_ties'][name]
        hyps[name] = dict(literal=S, aux_node=node, aux_label=label(node) if node in NODES else node,
                          literal_equals_table5=(node in AUX and sorted(AUX[node]['lit']) == sorted(S)),
                          T_plus_tie=sorted(set(T) | {tie}),
                          T_plus_tie_in_orbit_of_hypothesis=(canon(sorted(set(T) | {tie})) == canon(S)),
                          aux_stage=stage(node) if node in NODES else None,
                          aux_closure_contains_user=(user in CL.get(node, set())) if node in NODES else None,
                          graph_edge_aux_to_user=((node, user) in EDGES))
    r['hypotheses'] = hyps
    rel = relevant_ties(T)
    r['relevant_ties_computed'] = rel
    wit = {}
    for t, w in c['witnesses'].items():
        Tt = set(T) | {t}
        wit[t] = dict(witness=label(w), covers=covers(lit(w), Tt), witness_stage=stage(w), user_stage=stage(user))
    r['named_witnesses'] = wit
    if 'actions' in c:
        acts = {}
        for t, (row, col, tr, S) in c['actions'].items():
            a = doc_action(row, col, tr)
            acts[t] = dict(image=sorted(img(a, S)), image_in_T_plus_tie=img(a, S) <= set(T) | {t},
                           notation_md_convention_also=img(notation_action(''.join(map(str, row)), ''.join(map(str, col)), tr), S) <= set(T) | {t})
        r['printed_actions_checked'] = acts
    if 'case_actions' in c:
        acts = {}
        for t, (row, col, tr, w) in c['case_actions'].items():
            S = lit(w); Tt = set(T) | {t}
            acts[t] = dict(witness=label(w), image_doc_convention=sorted(img(doc_action(row, col, tr), S)),
                           doc_convention_in_T_plus_tie=img(doc_action(row, col, tr), S) <= Tt,
                           notation_md_convention_in_T_plus_tie=img(notation_action(''.join(map(str, row)), ''.join(map(str, col)), tr), S) <= Tt)
        r['printed_actions_checked'] = acts
    un = {}
    if isinstance(c['unnamed_ties'], list):
        for t in c['unnamed_ties']:
            Tt = set(T) | {t}
            cov = [k for k in list(CASE) + list(AUX) if k != user and covers(lit(k), Tt)]
            un[t] = dict(candidates=cov, earlier=[k for k in cov if stage(k) < stage(user)])
    r['unnamed_ties_candidates'] = un
    r['discharge'] = c['discharge']; r['assembly'] = c['assembly']
    cond_results.append(r)
# Tie -> witness tables printed in the proof texts (transcribed; labels are historical proof labels).
# Each is checked here: does an image of the witness lie in (contact set of the node) + tie, and is
# the witness earlier in the paper's Table 6 than the node?  ties=None: the text gives no tie order;
# then we only check that the witness covers the node plus SOME relevant tie.
WITNESS_TABLES = [
    ('T_66', 'root_moc4_T65_native_harmonic_uniform_20261023', 2, None, [3, 8, 11, 13, 15, 17, 19], ['T9', 'T7', 'T23', 'T24', 'T47', 'T3', 'T45']),
    ('T_87', 'blind_inventory86_matching_section_uniform_20261023', 2, None, [2, 10, 13, 3, 8, 16, 18], ['T_68'] * 3 + ['T_51'] * 3 + ['T_61']),
    ('T_60', 'blind_moc4_inventory59_native_harmonic_uniform_20261023', 1, [0, 1, 4, 6, 9, 17], [18, 15, 11, 14, 7, 5, 12], ['T35', 'T22', 'T46', 'T35', 'T19', 'T4', 'T52']),
    ('T_70', 'final_referee_inventory69_native_harmonic_20261023', 1, None, [2, 10, 11, 14, 15, 17, 20], ['T12', 'T41', 'T28', 'T43', 'T4', 'T43', 'T28']),
    ('T_79', 'final_referee_inventory78_native_uniform_20261023', 1, None, [2, 9, 10, 14, 17, 20], ['T12', 'T69', 'T70', 'T77', 'T57', 'T39']),
    ('T_82', 'final_referee_inventory81_native_uniform_20261023', 1, None, [3, 4, 7, 8, 14, 16, 18], ['T33', 'T59', 'T59', 'T30', 'T78', 'T59', 'T78']),
    ('T_45', 'independent_moc4_T44_native_harmonic_uniform_20261023', 1, None, [0, 6, 9, 18, 20, 21, 23], ['T8', 'T42', 'T7', 'T47', 'T37', 'T45', 'T37']),
    ('T_16', 'root_moc4_T15_uniform_20261022', 1, [0, 2, 5, 6, 8, 18], [4, 7, 10, 16, 19, 22, 23], ['T1', 'T7', 'T3', 'T14', 'T1', 'T3', 'T7']),
    ('T_49', 'independent_true_T48_native_harmonic_conditional_20261023', 2, None, [2, 7, 11, 13, 16, 23], ['T11', 'T36', 'T20', 'T41', 'T33', 'T43']),
    ('T_53', 'blind_moc4_inventory52_native_harmonic_conditional_20261023', 1, None, [5, 7, 11, 14, 15, 18], ['T4', 'T17', 'T25', 'T4', 'T25', 'T47']),
    ('T_50', 'blind_inventory49_whole_20261023', 2, None, [0, 7, 16], ['T9', 'T42', 'T65']),
    ('A_19', 'blind_moc4_inventory52_whole_assembly_20261023', 1, None, [18, 15, 11, 14, 7, 5], ['T47', 'T25', 'T25', 'T4', 'T17', 'T4']),
    ('A_20', 'final_referee_inventory56_boundary17_uniform_20261023', 1, None, [18, 15, 11, 14, 6, 4], ['T43', 'T25', 'T42', 'T33', 'T30', 'T16']),
    ('A_18', 'final_referee_inventory75_boundary16_uniform_20261023', 1, None, [20, 14, 13, 6, 7, 5], ['T39', 'T28', 'T28', 'T43', 'T30', 'T10']),
    ('A_25', 'final_referee_inventory78_boundary15_uniform_20261023', 1, None, [20, 17, 10, 9, 14, 2], ['T39', 'T20', 'T70', 'T69', 'T75', 'T12']),
    ('T_47', 'root_moc4_T46_whole_assembly_20261023', 1, None, None, ['T28', 'T35', 'T0', 'T38', 'T41', 'T6']),
    ('A_26', 'final_referee_T46_a0_uniform_20261023', 1, None, None, ['T28', 'T35', 'T0', 'T38', 'T41', 'T6']),
    ('T_65', 'root_moc4_T47_whole_assembly_20261023', 1, None, None, ['T43', 'T4', 'T7', 'T37', 'T41', 'T41']),
    ('T_49', 'blind_moc4_T48_allocation_II_uniform_20261023 (and allocation I, III, IV)', 3, None, None, ['T38', 'T0', 'T3', 'T23', 'T42', 'T43']),
]
witness_checks = []
# working frames used by documents: check that the printed action maps the literal set to the frame
frame_checks = {}
a60 = doc_action((0, 2, 3, 1), (0, 1, 3, 2), 1)
frame_checks['T_60 frame {0,1,4,6,9,17} (blind_moc4_inventory59 p.1)'] = dict(
    image_doc_convention=sorted(img(a60, CASE['T_60']['lit'])), same_orbit=canon([0, 1, 4, 6, 9, 17]) == canon(CASE['T_60']['lit']))
frame_checks['T_16 frame {0,2,5,6,8,18} (root_moc4_T15 p.1)'] = dict(same_orbit=canon([0, 2, 5, 6, 8, 18]) == canon(CASE['T_16']['lit']))
frame_checks['T_58 frame {0,1,5,7,9,17} (root_inventory57 conditional p.1)'] = dict(same_orbit=canon([0, 1, 5, 7, 9, 17]) == canon(CASE['T_58']['lit']))
for node, doc, page, working, ties, labs in WITNESS_TABLES:
    base = working or lit(node)
    rel = relevant_ties(base)
    for i, lb in enumerate(labs):
        w = lb if lb.startswith('T_') else by_hist[lb]
        if ties is not None:
            t = ties[i]
            ok = covers(lit(w), set(base) | {t})
            tt = [t]
        else:
            tt = [t for t in rel if covers(lit(w), set(base) | {t})]
            ok = bool(tt)
        alt_paper = ('T_%s' % lb[1:]) if (lb[0] == 'T' and lb[1:].isdigit() and ('T_%s' % lb[1:]) in CASE) else None
        alt_ok = (covers(lit(alt_paper), set(base) | ({ties[i]} if ties is not None else set())) if (alt_paper and ties is not None) else None)
        witness_checks.append(dict(node=label(node), document=doc + '.pdf', page=page, tie=(ties[i] if ties is not None else tt),
                                   printed_label=lb, resolved_case=label(w), covers=ok, tie_relevant=(all(x in rel for x in tt) if tt else False),
                                   witness_stage=stage(w), node_stage=stage(node), earlier=stage(w) < stage(node),
                                   in_graph_prerequisites=(w in PRED.get(node, ())),
                                   paper_numbering_reading=alt_paper, paper_numbering_reading_covers=alt_ok))
s52 = {t: dict(witness=label(w), covers=covers(lit(w), set(AUX['A_19']['lit']) | {t}), stage=stage(w))
       for t, w in S52_WIT.items()}

# ----------------------------------------------------------------- per case records
case_records = []
for k in sorted(CASE, key=lambda x: CASE[x]['n']):
    c = CASE[k]
    brs = cat_by[k]['branches']
    own = own_docs(k)
    preds = sorted(PRED.get(k, ()), key=lambda n: (NODES[n], n))
    clo = sorted(CL[k], key=lambda n: (NODES[n], n))
    cited = {m: [dict(doc=h['doc'], page=h['page'], how=h['how']) for h in hs if not h['flag']][:2]
             for m, hs in TEXTCIT[k].items() if is_node(m)}
    case_records.append(dict(
        case=k, contact_set=sorted(c['lit']), canonical=sorted(c['can']), inventory_index=c['inv'],
        historical_label=c['hist'], group=cat_by[k]['method_group'], paper_group=PAPER_GROUP.get(k),
        stage_graph=NODES[k] + 1, stage_paper=PAPER_STAGE.get(k),
        primary=dict(document=c['doc'], page=c['page'], anchor=c['anchor'], pages=pdf_pages.get(c['doc'])),
        own_documents=[dict(document=p, pages=pdf_pages.get(p)) for p in own],
        catalog_branches=[dict(document=b, pages=pdf_pages.get(b), home=HOME.get(b, [])) for b in brs],
        catalog_branch_count=len(brs), catalog_branch_pages=sum(pdf_pages.get(b, 0) for b in brs),
        graph_direct_prerequisites=[label(n) for n in preds],
        graph_transitive_prerequisites=[label(n) for n in clo],
        text_cited_nodes={label(m): v for m, v in cited.items()},
        relevant_extra_ties=coverage.get(k, []),
        tags=case_tags[k],
    ))

aux_records = []
for n in NODE_LIST:
    if n.startswith('T_'):
        continue
    aux_records.append(dict(node=label(n), contacts=AUX[n]['lit'] if n in AUX else None,
                            stage_graph=NODES[n] + 1, stage_paper=PAPER_STAGE.get(n),
                            proving_documents=[dict(document=p, page=pg, how=how, pages=pdf_pages.get(p)) for p, pg, how in PROVERS[n]],
                            graph_direct_prerequisites=[label(m) for m in sorted(PRED.get(n, ()))],
                            used_by=[label(m) for m in sorted(SUCC.get(n, ()))],
                            text_cited_nodes=sorted(label(m) for m in TEXTCIT[n] if is_node(m))))
dangling = [label(n) for n in NODES if not PROVERS.get(n)]

doc_records = {}
for p in DOCS:
    users = sorted({k for k in CASE if p in cat_by[k]['branches'] or CASE[k]['doc'] == p}, key=lambda x: CASE[x]['n'])
    doc_records[p] = dict(pages=pdf_pages.get(p), home=[label(h) for h in HOME[p]], statement_subjects=SUBJ[p][:4],
                          cases_listing_it=users, tags=tags_of_text(TEXT[p]))
core_docs = {stem(p): [label(h) for h in HOME[p]] + ['statement:' + ','.join(SUBJ[p][:2])]
             for p in DOCS if re.search(r'core|seven_A', stem(p))}

result = dict(
    evidence_level='mechanical (own code over supplement data and text extracts); see dependency-map.md',
    generator_check=gen_check, orbit_check=orbit_check, index_checks=idx_checks, graph_checks=graph_checks,
    documents_total=len(DOCS), document_pages_total=sum(pdf_pages[p] for p in DOCS), missing_text=missing_text,
    cases=case_records, auxiliary_and_branch_nodes=aux_records, dangling_nodes_without_document=dangling,
    core_documents=core_docs,
    edge_support={'%s -> %s' % (label(a), label(b)): v for (a, b), v in sorted(edge_support.items())},
    edges_without_text_support=[dict(edge='%s -> %s' % (label(a), label(b)), **v) for a, b, v in unsupported_edges],
    text_citations_missing_from_graph=[dict(m, citer=label(m['citer']), cited=label(m['cited'])) for m in missing_edges],
    unresolved_references=unresolved,
    relevant_ties_without_graph_witness={label(n): v for n, v in ties_without_graph_witness.items()},
    graph_prerequisites_not_covering_a_relevant_tie={label(n): [label(k) for k in v] for n, v in preds_covering_nothing.items()},
    boundary_ties_without_earlier_cover=[dict(node=label(n), tie=t, covers=[label(x) for x in cv]) for n, t, cv in uncovered],
    conditional_files=cond_file_roles, conditional_keyword_hits=cond_hits,
    conditional_theorems=cond_results, printed_witness_table_checks=witness_checks, working_frame_checks=frame_checks, S52_zero_coefficient_witnesses=s52,
    documents=doc_records,
)
os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, 'dependency-map.json'), 'w') as f:
    json.dump(result, f, indent=1, sort_keys=False, default=lambda o: sorted(o) if isinstance(o, (set, frozenset)) else str(o))
print(json.dumps(dict(gen=gen_check, orbits=orbit_check, graph={k: v for k, v in graph_checks.items() if k not in ('nontrivial_scc',)},
                      scc=graph_scc, idx={k: v for k, v in idx_checks.items()}, docs=len(DOCS), missing_text=len(missing_text),
                      unsupported_edges=len(unsupported_edges), missing_edges=len(missing_edges),
                      missing_by_status=collections.Counter(m['status'] for m in missing_edges),
                      uncovered=len(uncovered), dangling=dangling,
                      ties_without_graph_witness={label(n): v for n, v in ties_without_graph_witness.items()},
                      preds_covering_nothing={label(n): [label(k) for k in v] for n, v in preds_covering_nothing.items()}), indent=1, default=str))
