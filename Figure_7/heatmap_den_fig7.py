import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import string
    
fig, ax = plt.subplots(1,3,figsize=(9,4.4))

f1='/lnZ_nondegen/model_1/g_r_lnZ/E0d_e0_147_hc_a_0/'
f2='/lnZ_nondegen/model_1/g_r_lnZ/E0d_e0_91_hc_a_0/'

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f3='E0d_-0p034kT_mu_-15kT/gr.txt'
fx=''.join((f1,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
for i in range(2000):    
    o11[i]=(o1[i])/1000

ax1=plt.subplot(1, 3, 1)
ax1.text(-0.15, 1.05, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
plt.title("Fixed $l=147$")
plt.xlabel('Distance s (kbp)',fontsize=12)
plt.ylabel('$g^{(2)}$(s)',fontsize=12)
plt.xlim(0, 0.5)
plt.ylim(0,10)
plt.plot(o11[0:1000], o2[0:1000], c = 'black',label=r"$\mu=-15$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_-10kT/gr.txt'
fx=''.join((f1,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label=r"$\mu=-10$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_0kT/gr.txt'
fx=''.join((f1,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'b',linestyle='dashed',label=r"$\mu=0$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_10kT/gr.txt'
fx=''.join((f1,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'red',label=r"$\mu=10$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_15kT/gr.txt'
fx=''.join((f1,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label=r"$\mu=15$")

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
#plt.legend(loc="upper right",bbox_to_anchor=(1.03, 1.04),fontsize=10)
print('panel A...')

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f3='E0d_-0p034kT_mu_-15kT/gr.txt'
fx=''.join((f2,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
for i in range(2000):    
    o11[i]=(o1[i])/1000

ax1=plt.subplot(1, 3, 2)
ax1.text(-0.15, 1.05, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
plt.title("Fixed $l=91$")
plt.xlabel('Distance s (kbp)',fontsize=12)
plt.ylabel('$g^{(2)}$(s)',fontsize=12)
plt.xlim(0, 0.5)
plt.ylim(0,10)
plt.plot(o11[0:1000], o2[0:1000], c = 'black',label=r"$\mu=-15$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_-10kT/gr.txt'
fx=''.join((f2,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label=r"$\mu=-10$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_0kT/gr.txt'
fx=''.join((f2,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'b',linestyle='dashed',label=r"$\mu=0$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_10kT/gr.txt'
fx=''.join((f2,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'r',label=r"$\mu=10$")

o1=[]; o2=[];
f3='E0d_-0p034kT_mu_15kT/gr.txt'
fx=''.join((f2,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label=r"$\mu=15$")

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="lower right",bbox_to_anchor=(0.4, -0.7),borderpad=0.25,labelspacing=0.25,ncol=2,fontsize=10) #bbox_to_anchor=(1.03, 1.04)
print('panel B...')

o1=[]; o2=[]; o3=[]; o4=[]; o5=[]; 
o3x=np.zeros(20)
o11 = np.zeros(2000)
f3='mean_dg_vs_mu_fig7/d1_mu_E0d.txt'
fx=''.join((f2,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
    o4.append(float(row[3]))
    o5.append(float(row[4]))
f.close()
for i in range(20):    
    o3x[i]=1/((o3[i])*0.0001)

ax1=plt.subplot(1, 3, 3)
ax1.text(-0.15, 1.05, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Fixed-size 91")
plt.xlabel('NRL',fontsize=12)
plt.ylabel(r'1/$\rho$',fontsize=12)
plt.xlim(80, 200)
plt.ylim(80,200)
#plt.scatter(o4, o3x, s=8, c = 'grey',label='Fixed $l=91$')
plt.plot(o4, o3x,'o-',c = 'grey',fillstyle='none',label='Fixed $l=91$')

o1=[]; o2=[]; o3=[]; o4=[]; o5=[]; 
o3x=np.zeros(20)
o11 = np.zeros(2000)
f3='mean_dg_vs_mu_fig7/d1_mu_E0d.txt'
fx=''.join((f1,f3))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
    o4.append(float(row[3]))
    o5.append(float(row[4]))
f.close()
for i in range(20):    
    o3x[i]=1/((o3[i])*0.0001)

#plt.scatter(o4, o3x, s=8, c = 'black',label='Fixed $l=147$')
plt.plot(o4, o3x,'o-',c = 'grey',label='Fixed $l=147$')

o1=[]; o2=[]; o3=[]; o4=[]; o5=[]; 
o3x=np.zeros(20)
o11 = np.zeros(2000)
f = open('/lnZ_nondegen/mu_vs_E0d/mean_dg_vs_mu_fig7/e0_-0p034/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
    o4.append(float(row[3]))
    o5.append(float(row[4]))
f.close()
for i in range(20):    
    o3x[i]=1/((o3[i])*0.0001)

plt.plot(o4, o3x,'o-',c = 'red',label='model-1')
#plt.scatter(o4, o3x, s=8, c = 'red',label='model-1')


xx=np.zeros(21)
for i in range(21):    
    xx[i]=80+6*i
yy=xx
plt.plot(xx, yy,linestyle='dashed',c = 'grey',label=r'$1/\rho=$NRL')
print('panel C...')

plt.legend(loc="lower right",bbox_to_anchor=(1.0, -0.75),borderpad=0.25,labelspacing=0.25,ncol=1,fontsize=10) #bbox_to_anchor=(1.03, 1.04)

fig.tight_layout()
fig.subplots_adjust(wspace=0.32, hspace=0.25,bottom=0.4)
plt.show()


  
 




