#!/usr/bin/env python3
"""Independent, title-free alpha and beta crystal views; NumPy/Matplotlib only."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Circle

from crystal_geometry import (load_data, periodic_distances, orthographic, coordinate_frame,
                              COLOURS, draw_cell, beta_layer, beta_height)

HERE=Path(__file__).resolve().parent
OUT=HERE/'output'
BLUE=COLOURS['I']
ORANGE=COLOURS['II']


def figure(projection=False):
    fig,ax=plt.subplots(figsize=(89/25.4,89/25.4))
    if projection:
        fig.subplots_adjust(left=.17,right=.97,bottom=.19,top=.97)
    else:
        fig.subplots_adjust(left=.025,right=.975,bottom=.10,top=.975)
    return fig,ax


def save(fig,stem):
    assert not fig.texts and all(not ax.get_title() for ax in fig.axes)
    fig.canvas.draw()
    if len(fig.axes)>1:
        try:
            from audit_panel_alignment import require_matplotlib_panel_alignment
        except ImportError:
            w,h=fig.get_size_inches()*72
            boxes=np.array([[a.get_position().x0*w,a.get_position().y0*h,
                             a.get_position().x1*w,a.get_position().y1*h] for a in fig.axes])
            deviations=[np.ptp(boxes[:,2]-boxes[:,0]),np.ptp(boxes[:,3]-boxes[:,1]),
                        np.ptp(boxes[:,1]),np.ptp(boxes[:,3]),
                        np.ptp(boxes[1:,0]-boxes[:-1,2])]
            if max(deviations)>1.5: raise RuntimeError('View rectangles failed alignment')
            (OUT/(stem+'-alignment.json')).write_text(json.dumps(
                {'verdict':'PASS','tolerance_pt':1.5,'plot_rectangles_pt':boxes.tolist()},indent=2)+'\n')
        else:
            require_matplotlib_panel_alignment(fig,json_out=OUT/(stem+'-alignment.json'),
                                              tolerance_pt=1.5,gutter_tolerance_pt=1.5,
                                              require_panel_labels=False,strict=True)
    fig.savefig(OUT/(stem+'.png'),dpi=600,facecolor='white')
    fig.savefig(OUT/(stem+'.pdf'),facecolor='white')
    fig.savefig(OUT/(stem+'.svg'),facecolor='white')
    fig.savefig(OUT/(stem+'-preview.png'),dpi=300,facecolor='white')
    if len(fig.axes)==1:
        (OUT/(stem+'-alignment.json')).write_text(json.dumps(
            {'applicable':False,'reason':'one independent view; no comparable panels'},indent=2)+'\n')
    plt.close(fig)


def project(frac,item):
    return orthographic(np.atleast_2d(frac)*item['lengths'])[0]


def span_for(item):
    corners=project(np.array([[0,0,0],[1,1,1],[0,1,0],[1,0,1]]),item)
    return float(np.ptp(corners,axis=0).max())


def alpha_views(item,y):
    original=item['frac']
    shifted=(original-original[1])%1
    np.testing.assert_allclose(shifted[1],[0,0,0],atol=1e-12)
    np.testing.assert_allclose(shifted[0],[0,2*y,.5],atol=1e-12)
    np.testing.assert_allclose(periodic_distances(original,item['lengths']),
                               periodic_distances(shifted,item['lengths']),atol=1e-12,rtol=0)
    translated=dict(item,frac=shifted)
    fig,ax=figure()
    draw_cell(ax,translated,'alpha')
    span=span_for(translated)
    for index in (1,0):
        xy=project(shifted[index],translated)[0]
        ax.add_patch(Circle(xy,.041*span,facecolor='none',edgecolor=ORANGE,lw=1,zorder=8))
    origin=project(shifted[1],translated)[0]
    # A short leader identifies the chosen origin, without changing atom symbols.
    position=origin+np.array([.20,-.40])*span
    ax.text(*position,r'$(0,0,0)$',fontsize=8,color=ORANGE,ha='left',va='top')
    ax.plot([origin[0],position[0]],[origin[1],position[1]+.025*span],
            lw=.6,color=ORANGE,zorder=8)
    save(fig,'alpha-cell-origin')

    fig,ax=figure(projection=True)
    coordinate_frame(ax,r'$r_b/b$',r'$r_c/c$')
    # All four sites belong to the same 4c Wyckoff orbit: identical filled circles.
    for point in shifted:
        ax.plot(point[1],point[2],marker='o',ls='none',ms=5.5,
                markerfacecolor=BLUE,markeredgecolor=BLUE,markeredgewidth=.8)
    for index in (1,0):
        point=shifted[index]
        ax.plot(point[1],point[2],marker='o',ls='none',ms=7.2,
                markerfacecolor='none',markeredgecolor=ORANGE,markeredgewidth=1)
    two_y=2*y
    ax.plot([two_y,two_y],[0,.5],ls=':',lw=.7,color=ORANGE)
    for x in (0,two_y):
        ax.plot([x,x],[-.075,0],lw=.55,color=ORANGE)
    ax.annotate('',xy=(two_y,-.055),xytext=(0,-.055),arrowprops=dict(
        arrowstyle='<->',lw=.8,color=ORANGE,shrinkA=0,shrinkB=0))
    ax.text(two_y/2,.035,r'$2y$',color=ORANGE,fontsize=9,ha='center',va='bottom')
    save(fig,'alpha-bc-2y')
    return shifted


def local_legend(fig,sites):
    handles=[Line2D([],[],marker='o',ls='none',ms=4,color=COLOURS[site],label=site)
             for site in sites]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.005),
               ncol=len(sites),frameon=False,fontsize=8,handletextpad=.35,columnspacing=1)


def beta_views(item,reps):
    fig,ax=figure()
    draw_cell(ax,item,'beta')
    local_legend(fig,('I','II','III','IV','V'))
    save(fig,'beta-cell-sites')

    fig,ax=figure(projection=True)
    beta_layer(ax,item,reps)
    ax.set_xlabel(r'$r_a/a$',fontsize=8)
    ax.set_ylabel(r'$r_b/a$',fontsize=8)
    # This view has no III atoms; show only the classes actually present.
    local_legend(fig,('I','II','IV','V'))
    save(fig,'beta-ab-z0')

    fig,ax=figure(projection=True)
    beta_height(ax,item,reps)
    ax.set_xlabel(r'$r_a/a$',fontsize=8)
    ax.set_ylabel(r'$r_c/c$',fontsize=8)
    # Direct identification makes this projection independent of a remote legend.
    x,_,z=item['frac'][reps['III']]
    ax.annotate('III',xy=(x,z),xytext=(x+.16,z+.05),fontsize=8,color=COLOURS['III'],
                arrowprops=dict(arrowstyle='-',lw=.5,color=COLOURS['III']))
    save(fig,'beta-ac-III')


def beta_composite(item,reps):
    fig,axes=plt.subplots(1,3,figsize=(180/25.4,80/25.4))
    fig.subplots_adjust(left=.065,right=.975,bottom=.16,top=.96,wspace=.34)
    draw_cell(axes[0],item,'beta')
    beta_layer(axes[1],item,reps)
    beta_height(axes[2],item,reps)
    handles=[Line2D([],[],marker='o',ls='none',color=COLOURS[s],label=s,ms=4) for s in COLOURS]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.01),
               ncol=5,frameon=False,fontsize=8)
    save(fig,'beta-composite')


def main():
    global OUT
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True,help='local archived structure JSON; not bundled')
    parser.add_argument('--out',type=Path,default=HERE/'output')
    args=parser.parse_args()
    OUT=args.out
    OUT.mkdir(parents=True,exist_ok=True)
    source,structures,reps,params,checks=load_data(args.data)
    plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
                         'font.size':8,'pdf.fonttype':42,'svg.fonttype':'none','axes.linewidth':.6})
    shifted=alpha_views(structures['alpha'],params['alpha_y'])
    beta_views(structures['beta'],reps)
    beta_composite(structures['beta'],reps)
    fig,ax=figure()
    draw_cell(ax,structures['gamma'],'gamma')
    save(fig,'gamma-cell')
    metadata={'alpha_wyckoff':'4c; same filled-circle marker for every atom',
              'alpha_shift_origin':structures['alpha']['frac'][1].tolist(),
              'alpha_shifted_fractional':shifted.tolist(),'y':params['alpha_y'],
              'two_y':2*params['alpha_y'],'pair_distances_preserved':True,
              'orange_rings':'dimension endpoints, not separate Wyckoff classes',
              'beta_view_sites':{'cell':['I','II','III','IV','V'],
                                 'ab_z0':['I','II','IV','V'],'ac':['III']},
              'beta_view_unique_atom_counts':{'cell':30,'ab_z0':11,'ac_III':8},
              'source_checks':checks,'titles':False,'captions':False,
              'data_repository':source.get('repository'),
              'coordinate_source_commit':source.get('commit')}
    (OUT/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n',encoding='utf-8')
    print('Exported 2 alpha views, beta composite + 3 local-label views, and gamma cell.')


if __name__=='__main__': main()
