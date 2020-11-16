import numpy as np 
import matplotlib.pyplot as plt
import context
import settings

alpha=0.1

cmap_r=settings.colour_custom([194,24,91,1])
cmap_g=settings.colour_custom([76,175,80,1])
fig,axs=plt.subplots(1,2)
axs[0].set_title('a) BBNBI deposition')
axs[1].set_title('b) NUBEAM deposition')
axs[0].axvline(np.log10(alpha),color='r')
axs[1].axvline(np.log10(alpha),color='r')
axs[0].set_xlim(-5.,0.)
axs[1].set_xlim(-5.,0.)
bar_width=0.5
bar_pos=0.

ax=0
#****W03 KS tests****

'''

DFNs=['LOCUST GC', 'LOCUST GC truncated']
D=0.012493770445736324
P=3.54444360080373e-05
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width  
#reject? = True

DFNs=['LOCUST GC $\mathrm{ln}(\Lambda)_{\mathrm{ASCOT}}$ truncated', 'ASCOT GC']
D=0.01104420960511679
P=0.0003873883302903145
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True

'''
DFNs=['LOCUST GC $\mathrm{ln}(\Lambda)_{\mathrm{ASCOT}}$', 'ASCOT GC']
D=0.009653887401735117
P=0.0029114804077050036
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True
'''

DFNs=['LOCUST GC', 'ASCOT GC']
D=0.011639537530837019
P=0.00015034157070497854
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True

DFNs=['LOCUST FO $\mathrm{ln}(\Lambda)_{\mathrm{ASCOT}}$ truncated', 'ASCOT FO']
D=0.01319928945242077
P=9.948721417566468e-06
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True
'''

DFNs=['LOCUST FO $\mathrm{ln}(\Lambda)_{\mathrm{ASCOT}}$', 'ASCOT FO']
D=0.0007502884315115921
P=0.9999999999999996
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = False

DFNs=['LOCUST FO', 'ASCOT FO']
D=0.003340506079260286
P=0.8290390980567248
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = False

DFNs=['LOCUST FO', 'LOCUST FO truncated']
D=0.01333957306634509
P=7.663741002646893e-06
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True

DFNs=['LOCUST FO','LOCUST GC']
D=0.0025483302554601406
P=0.9767588466089407
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=2.*bar_width
#reject? = False

axs[ax].set_ylim(-0.3,bar_pos+0.7)
axs[ax].text(x=0.73,y=0.96,s='fail',fontsize=20,horizontalalignment='left',transform=axs[ax].transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.))
axs[ax].text(x=0.89,y=0.96,s='pass',fontsize=20,horizontalalignment='right',transform=axs[ax].transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.))


















bar_pos=0.
ax+=1
#****W04 KS tests****

'''
DFNs=['LOCUST GC', 'LOCUST GC truncated']
D=0.012433698824410033
P=3.9366759931675586e-05
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True

DFNs=['LOCUST GC $\mathrm{ln}(\Lambda)_{\mathrm{ASCOT}}$ truncated', 'ASCOT GC']
D=0.010631825131573111
P=0.000724842753918113
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True
'''

DFNs=['LOCUST GC $\mathrm{ln}(\Lambda)_{\mathrm{ASCOT}}$', 'ASCOT GC']
D=0.005828826690140287
P=0.1847087421907302
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True

'''

DFNs=['LOCUST GC', 'ASCOT GC']
D=0.008157192076601594
P=0.01886126541742162
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=3.*bar_width
#reject? = True
'''

DFNs=['LOCUST GC truncated', 'NUBEAM GC']
D=0.004287605913149584
P=0.5398293497378136
axs[ax].barh(bar_pos,np.log10(D),bar_width,color=cmap_g(0.))
axs[ax].barh(bar_pos+bar_width,np.log10(P),bar_width,color=cmap_r(0.))
label = axs[ax].annotate(
    ' - '.join([DFN for DFN in DFNs]), xy=(-1.1, bar_pos+2.*bar_width), xytext=(0, 0),
    textcoords="offset points",
    horizontalalignment='right', verticalalignment='center',
     clip_on=True, fontsize='x-small',color='black')
bar_pos+=2.*bar_width
#reject? = False

axs[ax].set_ylim(-0.3,bar_pos+0.3)

axs[ax].text(x=0.73,y=0.96,s='fail',fontsize=20,horizontalalignment='left',transform=axs[ax].transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.))
axs[ax].text(x=0.89,y=0.96,s='pass',fontsize=20,horizontalalignment='right',transform=axs[ax].transAxes,color=settings.colour_custom(rgba=[100,100,100,1])(0.))

fig.set_tight_layout(True)
axs[0].set_yticks([], [])
axs[1].set_yticks([], [])
axs[1].legend(['log alpha','log KS statistic','log P value'],loc='lower left')
plt.show()