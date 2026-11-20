#!/usr/bin/env python3
"""Validated conventional-cell coordinates and fractional projections."""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Rectangle

HERE = Path(__file__).resolve().parent
COLOURS = {'I': '#0072B2', 'II': '#D55E00', 'III': '#009E73',
           'IV': '#CC79A7', 'V': '#B17800'}
PHASES = ('alpha', 'beta', 'gamma')


def periodic_distances(frac, lengths):
    delta = frac[:, None, :] - frac[None, :, :]
    delta -= np.rint(delta)
    return np.linalg.norm(delta * lengths, axis=-1)


def load_data(path):
    raw = path if isinstance(path, dict) else json.loads(Path(path).read_text(encoding='utf-8'))
    structures, checks = {}, {}
    for phase in PHASES:
        item = raw['structures'][phase]
        f = np.array(item['fractional_raw'], dtype=float)
        lengths = np.array(item['cell_lengths_A'], dtype=float)
        expected = {'alpha': 4, 'beta': 30, 'gamma': 2}[phase]
        if f.shape != (expected, 3) or not np.isfinite(f).all():
            raise ValueError(f'{phase}: incomplete/nonfinite atom coordinates')
        if lengths.shape != (3,) or not np.all(lengths > 0):
            raise ValueError(f'{phase}: invalid orthogonal cell lengths')
        if not np.isclose(np.prod(lengths), item['cell_volume_A3'], rtol=5e-7):
            raise ValueError(f'{phase}: lengths disagree with archived volume')
        transformed = f.copy()
        if phase == 'beta':
            transformed[:, 1] *= -1
        transformed %= 1
        before, after = periodic_distances(f, lengths), periodic_distances(transformed, lengths)
        np.testing.assert_allclose(before, after, atol=1e-12, rtol=0)
        checks[phase] = {'atoms': expected, 'volume_consistent': True,
                         'max_periodic_distance_change_A': float(np.max(np.abs(before-after)))}
        structures[phase] = dict(item, frac=transformed, lengths=lengths)
    b = structures['beta']
    representatives = {'I': 0, 'II': 2, 'III': 26, 'IV': 21, 'V': 12}
    f = b['frac']
    # Check actual representatives, rather than attaching parameters to atom indices
    # whose coordinates happen to belong to a different member of the same orbit.
    ii, iii, iv, v = [f[representatives[s]] for s in ('II', 'III', 'IV', 'V')]
    np.testing.assert_allclose(ii, [ii[0], 1-ii[0], 0], atol=1e-10)
    np.testing.assert_allclose(iii[:2], [iii[0], iii[0]], atol=1e-10)
    np.testing.assert_allclose([iv[2], v[2]], [0, 0], atol=1e-10)
    a = structures['alpha']['frac']
    y = a[0, 1]
    ideal = np.array([[0,y,.25], [0,1-y,.75], [.5,.5+y,.25], [.5,.5-y,.75]])
    np.testing.assert_allclose(a, ideal, atol=2e-6)
    parameters = {'alpha_y': float(y), 'x_II': float(ii[0]), 'x_III': float(iii[0]),
                  'z_III': float(iii[2]), 'x_IV': float(iv[0]), 'y_IV': float(iv[1]),
                  'x_V': float(v[0]), 'y_V': float(v[1])}
    return raw, structures, representatives, parameters, checks



def orthographic(points):
    elev, azim = np.radians([20., 30.])
    depth = np.array([np.cos(elev)*np.cos(azim), np.cos(elev)*np.sin(azim), np.sin(elev)])
    right = np.array([-np.sin(azim), np.cos(azim), 0.])
    up = np.cross(depth, right)
    return points @ np.array([right, up]).T, points @ depth


def images(frac, labels):
    rows = []
    for f, label in zip(frac, labels):
        for shift in itertools.product((-1, 0, 1), repeat=3):
            point = f+shift
            if np.all(point >= -1e-8) and np.all(point <= 1+1e-8):
                rows.append((point, label))
    return rows



def coordinate_frame(ax, xlabel, ylabel):
    ax.set_aspect("equal", adjustable="box")
    # All diagrams use a unit fractional-coordinate square so plot areas align.
    ax.add_patch(Rectangle((0,0),1,1, fill=False, edgecolor='#A0A6AD', lw=.8))
    ax.set(xlim=(-.15,1.15), ylim=(-.15,1.15))
    ax.set_xticks([0,.5,1], ['0','1/2','1'])
    ax.set_yticks([0,.5,1], ['0','1/2','1'])
    ax.set_xlabel(xlabel, fontsize=8, labelpad=3)
    ax.set_ylabel(ylabel, fontsize=8, labelpad=3)
    ax.spines[['top','right']].set_visible(False)
    ax.spines[['left','bottom']].set_bounds(0,1)
    ax.tick_params(labelsize=8, length=2, pad=2)


def draw_cell(ax, item, phase):
    replicas = images(item['frac'], item['labels'])
    xy, depth = orthographic(np.array([r[0] for r in replicas])*item['lengths'])
    frac = np.array(list(itertools.product((0.,1.), repeat=3)))
    corners, _ = orthographic(frac*item['lengths'])
    span = float(np.ptp(corners,axis=0).max())
    for i,j in itertools.combinations(range(8),2):
        if np.sum(frac[i] != frac[j]) == 1:
            ax.plot(corners[[i,j],0], corners[[i,j],1],color='#9DA5AD',lw=.8,zorder=0)
    for rank,i in enumerate(np.argsort(depth)):
        colour = COLOURS[replicas[i][1]] if phase=='beta' else COLOURS['I']
        ax.add_patch(Circle(xy[i], .034*span, facecolor=colour,edgecolor='white',
                            lw=.45,zorder=2+rank*.01))
    centre = (corners.min(axis=0)+corners.max(axis=0))/2
    half = span*.73
    ax.set(xlim=(centre[0]-half,centre[0]+half),ylim=(centre[1]-half,centre[1]+half))
    ax.set_aspect('equal',adjustable='box');ax.axis('off')
    anchor=centre+np.array([-.58*span,-.57*span])
    for i,(name,offset) in enumerate(zip(('a','b','c'),((-.05,-.05),(.05,-.03),(0,.06)))):
        end=anchor+orthographic((np.eye(3)[i]*span*.19)[None,:])[0][0]
        ax.annotate('',xy=end,xytext=anchor,arrowprops=dict(arrowstyle='->',lw=.8,color='#444444'))
        ax.text(*(end+np.array(offset)*span),'$'+name+'$',fontsize=9,ha='center',va='center')


def beta_layer(ax,item,reps):
    coordinate_frame(ax,r'$x_{\rm frac}=r_a/a$',r'$y_{\rm frac}=r_b/a$')
    f=item['frac']
    for point,label in zip(f,item['labels']):
        if np.isclose(point[2],0,atol=1e-8):
            ax.plot(*point[:2],'o',ms=5.5,markerfacecolor=COLOURS[label],
                    markeredgecolor='white',markeredgewidth=.4)
    labels={'II':r'$x_{\rm II}$','IV':r'$(x_{\rm IV},y_{\rm IV})$',
            'V':r'$(x_{\rm V},y_{\rm V})$'}
    label_positions={'II':(.46,.69),'IV':(.50,.29),'V':(.48,.55)}
    for site in ('II','IV','V'):
        point=f[reps[site]];colour=COLOURS[site]
        ax.plot([0,point[0],point[0]],[point[1],point[1],0],color=colour,ls=':',lw=.65,zorder=0)
        position=label_positions[site]
        ax.text(*position,labels[site],va='center',ha='left',fontsize=8,color=colour)
        ax.plot([point[0],position[0]-.025],[point[1],position[1]],lw=.5,color=colour)


def beta_height(ax,item,reps):
    coordinate_frame(ax,r'$x_{\rm frac}=r_a/a$',r'$z_{\rm frac}=r_c/c$')
    for point,label in zip(item['frac'],item['labels']):
        if label=='III': ax.plot(point[0],point[2],'o',ms=5.5,color=COLOURS['III'])
    x,_,z=item['frac'][reps['III']]
    ax.plot([0,x,x],[z,z,0],color=COLOURS['III'],ls=':',lw=.7)
    ax.annotate('',xy=(x,.06),xytext=(0,.06),arrowprops=dict(
        arrowstyle='<->',lw=.8,color=COLOURS['III'],shrinkA=0,shrinkB=0))
    ax.text(x/2,.09,r'$x_{\rm III}a$',ha='center',va='bottom',fontsize=8,color=COLOURS['III'])
    ax.annotate('',xy=(-.06,z),xytext=(-.06,0),arrowprops=dict(
        arrowstyle='<->',lw=.8,color=COLOURS['III'],shrinkA=0,shrinkB=0))
    ax.annotate(r'$z_{\rm III}c$',xy=(-.06,z/2),xytext=(.03,.38),
                fontsize=8,color=COLOURS['III'],arrowprops=dict(arrowstyle='-',lw=.5,color=COLOURS['III']))
