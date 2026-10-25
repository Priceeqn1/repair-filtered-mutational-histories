#!/usr/bin/env python3
"""Generate a data-faithful vector overview (Paper I).

Requires: matplotlib (no input data file). All plotted points arise from the
finite weighted-composition recursion reproduced in this script.  The figure
retains the visual ideas of the originally generated concept artwork, but
replaces its illustrative event-count graph with normalized probabilities.
"""
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[1] / "figures"
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from repair_histories.model import event_distribution, two_hit_stats
R, G, EPS, K = 0.8, 1.3, 0.5, 10
RED, BLUE, INK, PALE = '#A83C45', '#276F9E', '#172532', '#EDF3F7'


s0, b0 = two_hit_stats(0)
s1, b1 = two_hit_stats(1)
assert abs(s0-5/21) < 1e-13
assert abs(s1-25/41) < 1e-13
assert abs(100*b1/b0-8.118811881188119) < .01
p0, p1 = np.array(event_distribution(K,0)), np.array(event_distribution(K,1))
assert abs(p0.sum()-1)<1e-13 and abs(p1.sum()-1)<1e-13
assert sum(np.arange(1,K+1)*p1) < sum(np.arange(1,K+1)*p0)

plt.rcParams.update({'font.family':'DejaVu Sans', 'mathtext.fontset':'dejavusans',
                     'pdf.fonttype':42, 'ps.fonttype':42, 'savefig.facecolor':'white'})
fig=plt.figure(figsize=(12.6, 5.25), facecolor='white')
# panel backgrounds are drawn in figure-normalized coords
shell=fig.add_axes([0,0,1,1]); shell.set_axis_off()
shell.set_xlim(0,1); shell.set_ylim(0,1)


def panel(x,w,label,title):
    y,h=.075,.805
    shell.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.006,rounding_size=.014',
                                  edgecolor='#A8B9C6',facecolor='white',linewidth=1.1))
    shell.add_patch(FancyBboxPatch((x+.007,.805),w-.014,.068,
                                  boxstyle='round,pad=.003,rounding_size=.013',
                                  edgecolor='none',facecolor=PALE))
    shell.text(x+.018,.836,label,fontsize=15,fontweight='bold',color=INK,va='center')
    shell.text(x+.065,.836,title,fontsize=11.7,fontweight='bold',color=INK,va='center')

shell.text(.5,.942,'Repair filters mutational histories—not only mutational burden',
           ha='center',va='center',fontsize=17.0,fontweight='bold',color=INK)
panel(.018,.297,'A','Same endpoint, two routes')
panel(.326,.338,'B','Two-hit quantitative effect')
panel(.675,.306,'C','General $k$-hit event ordering')

# --- Panel A: conceptual lineage schematic ---
def state(x,y,v,color):
    shell.add_patch(Circle((x,y),radius=.021,facecolor=color,edgecolor=INK,lw=1.15))
    shell.text(x,y,str(v),ha='center',va='center',color=INK,fontsize=12.3,fontweight='bold')

def arrow(a,b,color=INK,lw=2.3,style='-|>',rad=0):
    shell.add_patch(FancyArrowPatch(a,b,arrowstyle=style,mutation_scale=15,
                                   color=color,lw=lw,connectionstyle=f'arc3,rad={rad}'))

shell.text(.163,.737,'Direct double hit · 1 event',ha='center',fontsize=12.0,fontweight='bold',color=RED)
state(.084,.648,0,'#E9F0F4');state(.249,.648,2,'#F4D7DA')
arrow((.108,.648),(.222,.648),color=RED,lw=3.0)
shell.text(.163,.597,'Both alterations in one birth',ha='center',fontsize=10.2,color=INK)

shell.plot([.047,.283],[.550,.550],color='#D4DEE5',lw=1.0)
shell.text(.163,.504,'Sequential · 2 events',ha='center',fontsize=12.0,fontweight='bold',color=BLUE)
state(.070,.404,0,'#E9F0F4');state(.165,.404,1,'#D8EAF7');state(.260,.404,2,'#F4D7DA')
arrow((.092,.404),(.139,.404));arrow((.189,.404),(.233,.404))
shell.text(.115,.439,'birth 1',ha='center',fontsize=9.1,color=INK)
shell.text(.211,.439,'birth 2',ha='center',fontsize=9.1,color=INK)
arrow((.153,.380),(.086,.380),color=BLUE,lw=1.6,style='-|>',rad=-.58)
shell.text(.165,.307,'repair of the 1-hit state',ha='center',fontsize=10.8,color=BLUE,fontweight='bold')
shell.add_patch(FancyBboxPatch((.041,.135),.248,.107,boxstyle='round,pad=.008',
                              facecolor=PALE,edgecolor='none'))
shell.text(.165,.193,'Sequential accumulation has an',ha='center',fontsize=10.7,color=INK)
shell.text(.165,.162,'additional repairable intermediate state.',ha='center',fontsize=10.7,color=INK)

# --- Panel B: explicitly computed data with direct labels ---
xL,yTop=.358,.744
shell.text(.494,yTop,'Normalized two-hit abundance',ha='center',fontsize=11.5,color=INK,fontweight='bold')
def barchart(y, vals, ymax, labfmt):
    baseline=y
    shell.plot([.385,.609],[baseline,baseline],color=INK,lw=.9)
    for cx,v,c,nulab in zip((.439,.558),vals,(RED,BLUE),(r'$\nu=0$',r'$\nu=1$')):
        height=.123*v/ymax
        shell.add_patch(FancyBboxPatch((cx-.023,baseline),.046,height,boxstyle='square,pad=0',
                                      facecolor=c,edgecolor='none'))
        shell.text(cx,baseline+height+.015,labfmt(v),ha='center',va='bottom',fontsize=11.9,color=INK,fontweight='bold')
        shell.text(cx,baseline-.025,nulab,ha='center',va='top',fontsize=11.4,color=INK)

barchart(.555,(100.,100*b1/b0),105,lambda z:f'{z:.2f}%' if z!=100 else '100%')
shell.text(.494,.447,'Direct-route share among two-hit cells',ha='center',fontsize=11.5,fontweight='bold',color=INK)
barchart(.252,(100*s0,100*s1),75,lambda z:f'{z:.2f}%')
shell.add_patch(FancyBboxPatch((.361,.119),.264,.061,boxstyle='round,pad=.008',
                               facecolor=PALE,edgecolor='none'))
shell.text(.493,.151,r'$P_{\mathrm{direct}}(\nu)=\dfrac{\lambda_1+\nu}{\lambda_1+\nu+2\varepsilon r}$',
           ha='center',va='center',fontsize=15,color=INK)

# --- Panel C: exact conditional event count law (not illustrative proportions) ---
ax=fig.add_axes([.714,.280,.232,.433],facecolor='white')
x=np.arange(1,K+1)
# paired stems reduce occlusion; true normalized route distributions
ax.vlines(x-.13,0,p0,color=RED,lw=1.7,alpha=.9)
ax.vlines(x+.13,0,p1,color=BLUE,lw=1.7,alpha=.9)
ax.scatter(x-.13,p0,color=RED,s=26,zorder=3,label=r'Low repair ($\nu=0$)')
ax.scatter(x+.13,p1,color=BLUE,s=26,zorder=3,label=r'High repair ($\nu=1$)')
ax.set_xlim(.58,10.5); ax.set_ylim(0,.54)
ax.set_xticks(x);ax.set_yticks([0,.1,.2,.3,.4,.5]); ax.tick_params(labelsize=9,length=3,pad=2)
ax.set_ylabel(r'$\Pr(J=m\mid k)$',fontsize=11)
ax.set_xlabel('Mutation-bearing events, $m$',fontsize=10.4,labelpad=2)
ax.spines[['top','right']].set_visible(False)
ax.legend(loc='upper right',fontsize=8.2,frameon=False,handlelength=1.0,handletextpad=.4,
          borderaxespad=.1)
shell.text(.829,.751,r'Exact weighted routes: $k=10$',fontsize=11.2,
           ha='center',color=INK,fontweight='bold')
shell.add_patch(FancyBboxPatch((.690,.118),.277,.104,
                               boxstyle='round,pad=.008',facecolor=PALE,edgecolor='none'))
shell.text(.829,.190,'Higher repair shifts the full event-count',fontsize=10.6,ha='center',color=INK)
shell.text(.829,.162,'distribution toward fewer births.',fontsize=10.6,ha='center',color=INK)
shell.text(.829,.132,r'MLR ordering when $g+r\geq2$.',fontsize=10.1,ha='center',color=BLUE,fontweight='bold')

fig.text(.5,.030,
         r'Parameters: $r=0.8$, $g=1.3$, $\varepsilon=0.5$; $\nu\in\{0,1\}$.  All plotted values follow the rare-mutation route recursion.',
         ha='center',va='center',fontsize=10.8,color=INK)
for ext in ('pdf','png'):
    fig.savefig(ROOT/f'paper_I_three_panel_accurate.{ext}',dpi=220,transparent=False)
plt.close(fig)

# Save the exact probabilities used in panel C for independent inspection.
import csv
with (ROOT.parent / "data" / "figure_1_panel_C_event_counts.csv").open("w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["jump_count", "prob_nu0", "prob_nu1"])
    w.writerows(zip(range(1,K+1),p0,p1))
print('Pdirect:', 100*s0, 100*s1, 'burden ratio %:',100*b1/b0)
print('event-count sum:',p0.sum(),p1.sum())
print('k10 probabilities nu=0:',np.array2string(p0,precision=6))
print('k10 probabilities nu=1:',np.array2string(p1,precision=6))
