import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import string
import pandas as pd
import os
    
fig, ax = plt.subplots(2,2,figsize=(9,4.2))

o1=[]; o2=[]; o3=[]; o4=[];
Ns=57
amin=91
o11 = np.zeros(Ns)
for i in range(Ns):
    o11[i]=amin+i

f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/frag_91_147bp_gene.dat','r')
for row in f:
    row = row.split(' ')
    o2.append(float(row[0]))
f.close()

f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/frag_91_147bp_nong.dat','r')
for row in f:
    row = row.split(' ')
    o3.append(float(row[0]))
f.close()

i=0
sumo1=0
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/frag_91_147bp.dat','r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    sumo1=sumo1+o1[i];
f.close()
print(sumo1)

#f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/basic_nieb0_E0d_-3p08_mu_-10p5/frag_nie2_occ_0pxx.txt','r')
#for row in f:
#    row = row.split(' ')
#    o4.append(float(row[0]))
#f.close()

ax2=plt.subplot(1, 2, 1)
ax2.text(-0.18, 1.03, string.ascii_uppercase[0], transform=ax2.transAxes, size=20, weight='bold')
#plt.title("Non-degenerate")
plt.xlabel('Fragment length (bp)',fontsize=12)
plt.ylabel('Frequency',fontsize=13)
plt.axis([amin, 147, 0, 0.1])
#plt.yscale("log")


plt.plot(o11, o1, c = 'k')
plt.plot(o11, o2, c = 'm')
plt.plot(o11, o3, c = 'c')
#plt.plot(o11, o4, c = 'r')

plt.legend(('All', 'Gene', 'Non-gene')) #200U


#plt.xticks(np.arange(0, 400, 40))
#plt.yticks(np.arange(0, 1, 0.1))
plt.grid()
#plt.vlines(0, 0, 0.015, linestyles ="dashed", colors ="gray")
#plt.vlines(-0.05, 0, 0.015, linestyles ="dashed", colors ="red")
#plt.setp(ax2.get_yticklabels(), visible=False)

#ax2=plt.subplot(2, 2, 1)
#plt.setp(ax2.get_xticklabels(), visible=False)
#plt.setp(ax2.get_yticklabels(), visible=False)
#ax2.set_xticks([])
#ax2.set_yticks([])
#ax2=plt.subplot(2, 2, 3)
#plt.setp(ax2.get_xticklabels(), visible=False)
#plt.setp(ax2.get_yticklabels(), visible=False)
#ax2.set_xticks([])
#ax2.set_yticks([])


data3=np.load('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/pair_end/chrmap.npy')
s=data3.shape
print(s)
s=data3.sum()
print(s)
data4=data3[253:310,:]
o1_200u_sum=data4.sum()
data5=data4*(100/o1_200u_sum)
aa=data3.sum(axis=0)

ax4=plt.subplot(2, 2, 2)
ax4.text(-0.14, 1.06, string.ascii_uppercase[1], transform=ax4.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
#plt.xlabel('DNA (kbp)',fontsize=11)
plt.ylabel('Size (bp)',fontsize=12)
plt.xticks(np.arange(0, 2.01, 0.5))
plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(0, 2)
plt.ylim(91,147)
imx1=plt.imshow(data5,extent=[0,2,91,147], vmin=data5.min(), vmax=data5.max(), cmap='jet',aspect=0.015);
cax = ax4.inset_axes([1.03, 0.3, 0.02, 0.5])
cb=plt.colorbar(imx1,cax,shrink=0.4)
#plt.setp(ax4.get_xticklabels(), visible=False)
print(data5.max())
print(data5.min())

#data3=np.load('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/pair_end/chrmap.npy')

o1_200u=data5.sum(axis=0)
data5x=data5[0:7,:]
o2_200u=data5x.sum(axis=0)
data5x=data5[7:57,:]
o3_200u=data5x.sum(axis=0)

Ns=2000
o11 = np.zeros(2001)
for i in range(2001):
    o11[i]=i/1000
ax4=plt.subplot(2, 2, 4)
ax4.text(-0.14, 1.06, string.ascii_uppercase[2], transform=ax4.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Density',fontsize=12)

plt.xticks(np.arange(0, 2.01, 0.5))
plt.yticks(np.arange(0, 0.126, 0.025))
plt.ylim(0,0.125)
plt.xlim(0, 2)
#plt.yscale("log")
plt.plot(o11, o1_200u, c = 'k', label='All')
plt.plot(o11, o2_200u, c = 'b', label='Size$>140$')
plt.plot(o11, o3_200u, c = 'r', label=r'Size$\leq 140$')
plt.legend( loc="upper left",borderpad=0.25,labelspacing=0.25,fontsize=10)

fig.tight_layout()
plt.subplots_adjust(wspace=0.25, hspace=0.25)

plt.show()
 




