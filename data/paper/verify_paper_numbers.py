# -*- coding: utf-8 -*-
"""Verify every number in paper sections 4.3-4.6 against the data files in paper_data/.
Run from anywhere:  python verify_paper_numbers.py
Exit code 0 = all checks pass; 1 = at least one mismatch."""
import csv
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
FAIL = []
OK = 0

def check(label, computed, paper, tol=5.001e-4):
    global OK
    match = (computed == paper) if isinstance(paper, str) else (abs(computed - paper) <= tol)
    if match:
        OK += 1
        print(f'  [PASS] {label}: computed={computed!r}')
    else:
        FAIL.append(f'{label}: computed={computed!r} paper={paper!r}')
        print(f'  [FAIL] {label}: computed={computed!r} vs paper={paper!r}')

SEC44 = os.path.join(HERE, 'sec4_4_grid_convergence')

print('== 4.3 Table 3 (Re=100 row, g129, 40000 steps) ==')
j = json.load(open(os.path.join(SEC44, 'g129', 'time_convergence.json'), encoding='utf-8'))
c = [x for x in j['cases'] if x['steps'] == 40000][0]
check('u rel L2 %', round(c['u_relative_l2']*100, 3), 2.564)
check('u max err %', round(c['u_max_abs_percent'], 3), 1.984)
check('v rel L2 %', round(c['v_relative_l2']*100, 3), 3.982)
check('v max err %', round(c['v_max_abs_percent'], 3), 0.877)

print('== 4.4 Table 4 ==')
tc = {}
for g, N in (('g65', 65), ('g129', 129), ('g257', 257)):
    tc[N] = json.load(open(os.path.join(SEC44, g, 'time_convergence.json'), encoding='utf-8'))
c65 = [x for x in tc[65]['cases'] if x['steps'] == 40000][0]
c129 = c
c257 = [x for x in tc[257]['cases'] if x['steps'] == 40000][0]
for name, c_, vals in (('65x65', c65, (2.590, 2.201, 2.243, 0.483)),
                       ('129x129', c129, (2.564, 1.984, 3.982, 0.877)),
                       ('257x257', c257, (1.983, 1.552, 3.621, 0.780))):
    check(f'{name} u rel L2', round(c_['u_relative_l2']*100, 3), vals[0])
    check(f'{name} u max', round(c_['u_max_abs_percent'], 3), vals[1])
    check(f'{name} v rel L2', round(c_['v_relative_l2']*100, 3), vals[2])
    check(f'{name} v max', round(c_['v_max_abs_percent'], 3), vals[3])
check('order u 65->129', round(math.log2(c65['u_max_abs_percent']/c129['u_max_abs_percent']), 2), 0.15)
check('order u 129->257', round(math.log2(c129['u_max_abs_percent']/c257['u_max_abs_percent']), 2), 0.35)
check('order v 65->129', round(math.log2(c65['v_max_abs_percent']/c129['v_max_abs_percent']), 2), -0.86)
check('order v 129->257', round(math.log2(c129['v_max_abs_percent']/c257['v_max_abs_percent']), 2), 0.17)
check('div norm g65', round(c65['div_norm'], 4), 0.0097)
check('div norm g129', round(c129['div_norm'], 4), 0.0087)
check('div norm g257', round(c257['div_norm'], 4), 0.0063)
check('u_min g65', round(c65['u_min'], 3), -0.197)
check('u_min g129', round(c129['u_min'], 3), -0.193)
check('u_min g257', round(c257['u_min'], 3), -0.196)
s = {x['steps']: x for x in tc[257]['cases']}
check('g257 u relL2 @5k %', round(s[5000]['u_relative_l2']*100, 1), 13.6)
check('g257 u relL2 @10k %', round(s[10000]['u_relative_l2']*100, 1), 7.1)
check('g257 u relL2 @20k %', round(s[20000]['u_relative_l2']*100, 1), 2.4)
check('g257 u relL2 @40k %', round(s[40000]['u_relative_l2']*100, 1), 2.0)
check('g257 u_change 20k->40k %', round(s[40000]['u_change']*100, 1), 4.2)
check('g257 v_change 20k->40k %', round(s[40000]['v_change']*100, 1), 8.0)
check('g65 v_change %', round(c65['v_change']*100, 2), 0.70)
check('g129 u_change %', round(c129['u_change']*100, 2), 0.33)
check('g257 physical end time', 40000*0.0005, 20, tol=1e-9)
check('diffusion number g257 dt=0.001', round(0.02*0.001/(2/256)**2, 2), 0.33)
check('diffusion number g257 dt=0.0005', round(0.02*0.0005/(2/256)**2, 2), 0.16)

print('== 4.5 Table 5 (recomputed from field CSV) ==')
u = [[0.0]*129 for _ in range(129)]
v = [[0.0]*129 for _ in range(129)]
with open(os.path.join(HERE, 'sec4_5_flow_structure', 'field_re100_g129.csv'), encoding='utf-8') as f:
    for r in csv.DictReader(l for l in f if not l.startswith('#')):
        u[int(r['i'])][int(r['j'])] = float(r['u'])
        v[int(r['i'])][int(r['j'])] = float(r['v'])
dx = 2.0/128
psi = [[0.0]*129 for _ in range(129)]
for j in range(129):
    acc = 0.0
    for i in range(1, 129):
        acc += 0.5*(u[i-1][j] + u[i][j])*dx
        psi[i][j] = acc
best = (1e18, 0, 0)
for i in range(1, 128):
    for j in range(1, 128):
        if psi[i][j] < best[0]:
            best = (psi[i][j], i, j)
check('primary x*', round(best[2]*dx/2, 4), 0.6094)
check('primary y*', round(best[1]*dx/2, 4), 0.7266)
check('primary psi*', round(best[0]/2, 4), -0.0931)
umin = min((u[i][64], i) for i in range(129))
vmin = min((v[64][j], j) for j in range(129))
check('u_min', round(umin[0], 4), -0.1931)
check('v_min', round(vmin[0], 4), -0.2370)

def local_max(i0, i1, j0, j1):
    b = (-1e18, -1, -1)
    for i in range(i0, i1+1):
        for j in range(j0, j1+1):
            if psi[i][j] > b[0]:
                b = (psi[i][j], i, j)
    return b

_, bi, bj = local_max(1, 20, 1, 20)
check('BL1 x*', round(bj*dx/2, 4), 0.0313)
check('BL1 y*', round(bi*dx/2, 4), 0.0391)
bl1_psi = psi[bi][bj]/2
check('BL1 psi* (text 1.6e-6)', float(f'{bl1_psi:.1e}'), 1.6e-6, tol=0.051e-6)
_, bi, bj = local_max(1, 20, 108, 127)
check('BR1 x*', round(bj*dx/2, 4), 0.9453)
check('BR1 y*', round(bi*dx/2, 4), 0.0547)
br1_psi = psi[bi][bj]/2
check('BR1 psi* (text 9.1e-6)', float(f'{br1_psi:.1e}'), 9.1e-6, tol=0.051e-6)
check('dev x*', round(abs(0.60938-0.61621), 4), 0.0068)
check('dev y*', round(abs(0.72656-0.73730), 4), 0.0107)
check('dev psi %', round((0.1035212-0.0931236)/0.1035212*100, 1), 10.0)
check('dev u_min %', round((0.2140417-0.1930762)/0.2140417*100, 1), 9.8)
check('dev v_min %', round((0.253804-0.2370052)/0.253804*100, 1), 6.6)

print('== 4.6 Table 6 (recomputed from bous CSVs; u*=20u, v*=20v) ==')
bench = {3: (1.118, 3.649, 3.697), 4: (2.238, 16.178, 19.62), 5: (4.519, 34.81, 68.22)}
paper = {3: (1.166, 3.356, 4.005), 4: (2.633, 16.08, 23.18), 5: (6.397, 39.50, 93.31)}

def sf4(x):
    return 0.0 if x == 0 else round(x, -int(math.floor(math.log10(abs(x)))) + 3)

for ra, fn in ((3, 'bous_ra1e3.csv'), (4, 'bous_ra1e4.csv'), (5, 'bous_ra1e5.csv')):
    rows = [r for r in csv.DictReader(open(os.path.join(HERE, 'sec4_6_natural_convection', fn),
                                           encoding='utf-8'))]
    last, prev = rows[-1], rows[-2]
    nu = float(last['nusselt'])
    ust = 20*float(last['u_max_vc'])
    vst = 20*float(last['v_max_hc'])
    check(f'Ra=1e{ra} Nu', sf4(nu), paper[ra][0])
    check(f'Ra=1e{ra} u*', sf4(ust), paper[ra][1])
    check(f'Ra=1e{ra} v*', sf4(vst), paper[ra][2])
    rel = abs(nu - float(prev['nusselt']))/abs(nu)
    check(f'Ra=1e{ra} plateau rel change <1e-4 (got {rel:.2e})', rel < 1e-4, True)
    check(f'Ra=1e{ra} Nu error vs bench', round((nu-bench[ra][0])/bench[ra][0]*100, 1),
          {3: 4.3, 4: 17.7, 5: 41.6}[ra])
    check(f'Ra=1e{ra} u* err vs bench', round((ust-bench[ra][1])/bench[ra][1]*100, 1),
          {3: -8.0, 4: -0.6, 5: 13.5}[ra])
for beta, target in ((0.8875, 1e3), (8.875, 1e4), (88.75, 1e5)):
    check(f'Ra mapping beta={beta}', round(beta*8/(0.071*0.1)), round(target))

print()
print(f'RESULT: {OK} checks passed, {len(FAIL)} failed')
for f_ in FAIL:
    print('  FAIL:', f_)
sys.exit(1 if FAIL else 0)
