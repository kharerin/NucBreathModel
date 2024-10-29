import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import string
import pandas as pd
import os
    
#fig, ax1 = plt.subplots(5,2,figsize=(6.8,9))

fig = plt.figure(figsize=(8,12))
gs = gridspec.GridSpec(12, 4)
gs.update(wspace = 1.0, hspace = 0, bottom=0.055, top=0.95, left=0.09, right=0.912)

o1=[]; o2=[]; o3=[]; o4=[]; o5=[];
Ns=57
amin=91
o11 = np.zeros(Ns)
for i in range(Ns):
    o11[i]=amin+i

i=0
sumo1=0
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/frag_91_147bp.dat','r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    sumo1=sumo1+o1[i];
f.close()
print(sumo1)

f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/E0d_-3p6kT_mu_-9p2kT_amin91/frag_nie2_occ_0pxx.txt','r')
for row in f:
    row = row.split(' ')
    o3.append(float(row[0]))
f.close()

f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/basic_nieb0_E0d_-3p2_mu_-12p8/frag_nie2_occ_0pxx.txt','r')
for row in f:
    row = row.split(' ')
    o4.append(float(row[0]))
f.close()

f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/E0d_-5p4kT_mu_-11p6kT/frag_nie2_occ_0pxx.txt','r')
for row in f:
    row = row.split(' ')
    o5.append(float(row[0]))
f.close()
ax1 = plt.subplot(gs[0:2, 0:2])
ax1.text(-0.22, 1.1, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Non-degenerate")
plt.xlabel('Fragment length (bp)',fontsize=11)
plt.ylabel('Frequency',fontsize=11)
plt.axis([amin, 147, 0, 0.08])
#plt.yscale("log")


plt.plot(o11, o1, c = 'k',label='Expt.data')
plt.plot(o11, o4, c = 'r',label='model-1')
plt.plot(o11, o3, c = 'b',linestyle='-',label='model-2')
plt.plot(o11, o5[90:147], c = 'limegreen',linestyle='-',label='model-3')

#plt.legend(('Expt.data', '$l_{min}=91$,$l_{max}=147$,$l_k=0$', '$l_{min}=91$,$l_{max}=147$,$l_k=20$','$l_{min}=1$,$l_{max}=147$,$l_k=20$')) #200U
plt.legend(loc="upper left",borderpad=0.05,labelspacing=0.25,ncol=1,facecolor='w',framealpha=1,frameon=False,fontsize=10) #200U

#plt.xticks(np.arange(0, 400, 40))
#plt.yticks(np.arange(0, 1, 0.1))
#plt.grid()
#plt.vlines(0, 0, 0.015, linestyles ="dashed", colors ="gray")
#plt.vlines(-0.05, 0, 0.015, linestyles ="dashed", colors ="red")
#plt.setp(ax2.get_yticklabels(), visible=False)

data3=np.load('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/yeastsc_lnZ/H4_Input_MNase_200U_Rep1/pair_end/chrmap.npy')
s=data3.shape
print(s)
s=data3.sum()
print(s)
aa=data3.sum(axis=0)

data4=data3[253:310,:]
o1_200u_sum=data4.sum()
data5=data4*(100/o1_200u_sum)
data5x=data4*(100/o1_200u_sum)
s=data5.shape
print("sss",s)
o1_200u=data5.sum(axis=0)
data4=data3[253:260,:]
o2_200u_sum=data4.sum()
o2_200u=data4.sum(axis=0)
data4=data3[259:310,:]
o3_200u_sum=data4.sum()
o3_200u=data4.sum(axis=0)

Ns=2000
ax2=plt.subplot(gs[0:2, 2:4])
ax2.text(-0.22, 1.1, string.ascii_uppercase[1], transform=ax2.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
plt.xlabel('DNA (kbp)',fontsize=11)
plt.ylabel('Density',fontsize=11)

plt.xticks(np.arange(0, 1.1, 0.2))
plt.yticks(np.arange(0, 0.31, 0.1))
plt.ylim(0,0.3)
plt.xlim(0, 1)
#plt.yscale("log")
o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
plt.plot(o11, o1_200u[1000:2000], c = 'k',label='Expt.data') #plt.plot(o11[1000:2000], o1_200u[1000:2000], c = 'k')
#plt.plot(o11, o2_200u, c = 'b')
#plt.plot(o11, o3_200u, c = 'r')
#plt.legend(('All', 'Size>140', 'Size<141'))


d1=[]
d2=[]
d3=[]
data1 = np.zeros( (57,10000) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/basic_nieb0_E0d_-3p2_mu_-12p8/dyamapfull_h10_sigm8_nie2.txt','r')
for row in f:
    row = row.split(' ')
    d1.append(float(row[0]))
    d2.append(float(row[1]))
    d3.append(float(row[2]))
f.close()
k=0   
for i in range(2000):
    for j in range(57):
        data1[56-j][i]=d3[j+k]
    k=k+57
data1x=data1
sx=data1.sum()
data5=data1*(100/sx)
data5a=data5
aa=data5.sum(axis=0)
o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
plt.plot(o11, aa[25:1025], c = 'r',label='model-1') #plt.plot(o11[975:1975], aa[0:1000], c = 'r')


d1=[]
d2=[]
d3=[]
data1 = np.zeros( (57,10000) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/E0d_-3p6kT_mu_-9p2kT_amin91/dyamapfull_h10_sigm8_nie2.txt','r')
for row in f:
    row = row.split(' ')
    d1.append(float(row[0]))
    d2.append(float(row[1]))
    d3.append(float(row[2]))
f.close()
k=0   
for i in range(2000):
    for j in range(57):
        data1[56-j][i]=d3[j+k]
    k=k+57
data2x=data1
sx=data1.sum()
data5=data1*(100/sx)
data5b=data5
aa=data5.sum(axis=0)
o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
plt.plot(o11, aa[25:1025], c = 'b',linestyle='-',label='model-2') #plt.plot(o11[975:1975], aa[0:1000], c = 'r')

d1=[]
d2=[]
d3=[]
data1 = np.zeros( (147,10000) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/E0d_-5p4kT_mu_-11p6kT/dyamapfull_h10_sigm8_nie2.txt','r')
for row in f:
    row = row.split(' ')
    d1.append(float(row[0]))
    d2.append(float(row[1]))
    d3.append(float(row[2]))
f.close()
k=0   
for i in range(2000):
    for j in range(147):
        data1[146-j][i]=d3[j+k]
    k=k+147
data3x=data1[0:57,:]
sx=data1[0:57,:].sum()
data5=data1[0:57,:]*(100/sx)
data5c=data5
aa=data5.sum(axis=0)
o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
plt.plot(o11, aa[25:1025], c = 'limegreen',linestyle='-',label='model-3') #plt.plot(o11[975:1975], aa[0:1000], c = 'r')

d1=[]
d2=[]
d3=[]
data1 = np.zeros(10000 )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/g_r_lnZ/E0d_e0_147_hc_a_0/E0d_-0p41kT_mu_-54p8kT/dyamapfull_h10_sigm8_nie2.txt','r')
for row in f:
    row = row.split(' ')
    d1.append(float(row[0]))
    d2.append(float(row[1]))
    d3.append(float(row[2]))
f.close()
k=0   
for i in range(2000):
        data1[i]=d3[i]
sx=data1.sum()
aa=data1*(100/sx)
o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
plt.plot(o11, aa[25:1025], c = 'grey',linestyle='-',label='Fixed-147') #plt.plot(o11[975:1975], aa[0:1000], c = 'r')

plt.legend(loc="upper right",borderpad=0.05,labelspacing=0.25,ncol=2,facecolor='w',framealpha=1,frameon=False,fontsize=10) #200U

#------------------------------------------------------------------------------------------

#data4x=data5x[253:310,1000:2000]
#s=data4x.shape
#print(s)
#s=data4x.sum()
data4=data5x[:,1000:2000]
fs0=data4.sum()
fs1=data4[:,147:1000].sum()
fs2=data4[:,0:147].sum()
fs3=fs0-fs1
print(fs0,fs1,fs2,fs3)
fexpt=data4[:,0:147].sum(axis=1)
fexptx=np.zeros(57)
for i in range(57):
    fexptx[i]=fexpt[56-i]
print(s)
s=data4.sum()
print(s)

ax3=plt.subplot(gs[3:5, 0:2])
ax3.text(-0.22, 1.1, string.ascii_uppercase[2], transform=ax3.transAxes, size=20, weight='bold')
plt.title("Expt.data",fontsize=12)
#plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Size (bp)',fontsize=12)
plt.xticks(np.arange(0, 1.01, 0.2))
plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(0, 1)
plt.ylim(91,147)
imx1=plt.imshow(data4,extent=[0,1,91,147], vmin=0, vmax=0.015, cmap='jet',aspect=0.008);
#cax = ax4.inset_axes([1.03, 0.3, 0.02, 0.5])
#cb=plt.colorbar(imx1,cax,shrink=0.4)
plt.setp(ax3.get_xticklabels(), visible=False)
print(data4.max())
print(data4.min())

o1_200u=data4.sum(axis=0)
data5=data4[0:7,:]
o2_200u=data5.sum(axis=0)
data5=data4[7:57,:]
o3_200u=data5.sum(axis=0)

o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
ax5=plt.subplot(gs[5:7, 0:2])
#ax5.text(-0.12, 1.05, string.ascii_uppercase[2], transform=ax5.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Density',fontsize=12)

plt.xticks(np.arange(0, 1.1, 0.2))
plt.yticks(np.arange(0, 0.31, 0.1))
plt.ylim(0,0.3)
plt.xlim(0, 1)
#plt.yscale("log")
plt.plot(o11, o1_200u, c = 'k', label='All')
plt.plot(o11, o2_200u, c = 'b', label='Size$>140$')
plt.plot(o11, o3_200u, c = 'r', label='Size$\\leq 140$')
plt.legend( loc="upper right",borderpad=0.25,labelspacing=0.25,fontsize=10)
#Ns=57
#o11 = np.zeros(Ns)
#for i in range(Ns):
#    o11[i]=91+i
#axins = ax5.inset_axes([0.18, 0.56, 0.36, 0.40])
#axins.plot(o11, o1, c = 'k',label='Expt.data')
#axins.plot(o11, fexptx/fs2, c = 'r',label='model-1')

data4x=data5a[:,25:1025]
s=data4x.shape
print(s)
s=data4x.sum()
data4=data5a[:,25:1025]
print(s)
s=data4.sum()
print(s)

ax4=plt.subplot(gs[3:5, 2:4])
ax4.text(-0.22, 1.1, string.ascii_uppercase[3], transform=ax4.transAxes, size=20, weight='bold')
plt.title("model-1",fontsize=12)
#plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Size (bp)',fontsize=12)
plt.xticks(np.arange(0, 1.01, 0.2))
plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(0, 1)
plt.ylim(91,147)
imx1=plt.imshow(data4,extent=[0,1,91,147], vmin=0, vmax=0.015, cmap='jet',aspect=0.008);
#cax = ax4.inset_axes([1.03, 0.3, 0.02, 0.5])
#cb=plt.colorbar(imx1,cax,shrink=0.4)
plt.setp(ax4.get_xticklabels(), visible=False)
print(data4.max())
print(data4.min())

o1_200u=data4.sum(axis=0)
data5=data4[0:7,:]
o2_200u=data5.sum(axis=0)
data5=data4[7:57,:]
o3_200u=data5.sum(axis=0)

o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
ax6=plt.subplot(gs[5:7, 2:4])
#ax6.text(-0.25, 1.1, string.ascii_uppercase[2], transform=ax6.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Density',fontsize=12)

plt.xticks(np.arange(0, 1.1, 0.2))
plt.yticks(np.arange(0, 0.31, 0.1))
plt.ylim(0,0.3)
plt.xlim(0, 1)
#plt.yscale("log")
plt.plot(o11, o1_200u, c = 'k', label='All')
plt.plot(o11, o2_200u, c = 'b', label='Size$>140$')
plt.plot(o11, o3_200u, c = 'r', label='Size$\\leq 140$')
plt.legend( loc="upper right",borderpad=0.25,labelspacing=0.25,fontsize=10)


data4x=data5b[:,25:1025]
s=data4x.shape
print(s)
s=data4x.sum()
data4=data5b[:,25:1025]
fs0=data4.sum()
fs1=data4[:,147:1000].sum()
fs2=data4[:,0:147].sum()
fs3=fs0-fs1
print(fs0,fs1,fs2,fs3)
fexpt=data4[:,0:147].sum(axis=1)
fexptx=np.zeros(57)
for i in range(57):
    fexptx[i]=fexpt[56-i]
print(s)
s=data4.sum()
print(s)

ax7=plt.subplot(gs[8:10, 0:2])
ax7.text(-0.22, 1.1, string.ascii_uppercase[4], transform=ax7.transAxes, size=20, weight='bold')
plt.title("model-2",fontsize=12)
#plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Size (bp)',fontsize=12)
plt.xticks(np.arange(0, 1.01, 0.2))
plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(0, 1)
plt.ylim(91,147)
imx1=plt.imshow(data4,extent=[0,1,91,147], vmin=0, vmax=0.015, cmap='jet',aspect=0.008);
#cax = ax4.inset_axes([1.03, 0.3, 0.02, 0.5])
#cb=plt.colorbar(imx1,cax,shrink=0.4)
plt.setp(ax7.get_xticklabels(), visible=False)
print(data4.max())
print(data4.min())

o1_200u=data4.sum(axis=0)
data5=data4[0:7,:]
o2_200u=data5.sum(axis=0)
data5=data4[7:57,:]
o3_200u=data5.sum(axis=0)

o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
ax9=plt.subplot(gs[10:12, 0:2])
#ax9.text(-0.12, 1.05, string.ascii_uppercase[2], transform=ax9.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Density',fontsize=12)

plt.xticks(np.arange(0, 1.1, 0.2))
plt.yticks(np.arange(0, 0.31, 0.1))
plt.ylim(0,0.3)
plt.xlim(0, 1)
#plt.yscale("log")
plt.plot(o11, o1_200u, c = 'k', label='All')
plt.plot(o11, o2_200u, c = 'b', label='Size$>140$')
plt.plot(o11, o3_200u, c = 'r', label='Size$\\leq 140$')
plt.legend( loc="upper right",borderpad=0.25,labelspacing=0.25,fontsize=10)
#Ns=57
#o11 = np.zeros(Ns)
#for i in range(Ns):
#    o11[i]=91+i
#axins = ax9.inset_axes([0.18, 0.56, 0.36, 0.40])
#axins.plot(o11, o3, c = 'k',label='bulk')
#axins.plot(o11, fexptx/fs2, c = 'grey',linestyle='dashed',label='+1')
#x1, x2, y1, y2 = 90, 147, 0, 0.05
#axins.set_xlim(x1, x2) 
#axins.set_ylim(y1, y2)


data4x=data5c[:,25:1025]
s=data4x.shape
print(s)
s=data4x.sum()
data4=data5c[:,25:1025]
print(s)
s=data4.sum()
print(s)

ax8=plt.subplot(gs[8:10, 2:4])
ax8.text(-0.22, 1.1, string.ascii_uppercase[5], transform=ax8.transAxes, size=20, weight='bold')
plt.title("model-3",fontsize=12)
#plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Size (bp)',fontsize=12)
plt.xticks(np.arange(0, 1.01, 0.2))
plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(0, 1)
plt.ylim(91,147)
imx1=plt.imshow(data4,extent=[0,1,91,147], vmin=0, vmax=0.015, cmap='jet',aspect=0.008);
cax = ax1.inset_axes([2.4, -4.2, 0.025, 0.5])
cb=plt.colorbar(imx1,cax,shrink=0.4)
plt.setp(ax8.get_xticklabels(), visible=False)
print(data4.max())
print(data4.min())

o1_200u=data4.sum(axis=0)
data5=data4[0:7,:]
o2_200u=data5.sum(axis=0)
data5=data4[7:57,:]
o3_200u=data5.sum(axis=0)

o11 = np.zeros(1000)
for i in range(1000):
    o11[i]=i/1000
ax10=plt.subplot(gs[10:12, 2:4])
#ax10.text(-0.12, 1.05, string.ascii_uppercase[2], transform=ax10.transAxes, size=20, weight='bold')
#plt.title("Flat7",fontsize=11)
plt.xlabel('DNA (kbp)',fontsize=12)
plt.ylabel('Density',fontsize=12)

plt.xticks(np.arange(0, 1.1, 0.2))
plt.yticks(np.arange(0, 0.31, 0.1))
plt.ylim(0,0.3)
plt.xlim(0, 1)
#plt.yscale("log")
plt.plot(o11, o1_200u, c = 'k', label='All')
plt.plot(o11, o2_200u, c = 'b', label='Size$>140$')
plt.plot(o11, o3_200u, c = 'r', label='Size$\\leq 140$')
plt.legend( loc="upper right",borderpad=0.25,labelspacing=0.25,fontsize=10)

#plt.subplots_adjust(wspace=0, hspace=0)
#fig.tight_layout()

plt.show()
 




