_import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import math
import numpy as np 
import string

fig, ax = plt.subplots(3,3,figsize=(9,8))

mu=[]
E0=[]
O=[]
pO=[]
lm=[]
data1 = np.zeros( (6, 575) )
data4 = np.zeros( (6, 575) )
data5 = np.zeros( (6, 575) )
f = open('/lnZ_nondegen/model_1/density_amin91_vs_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    O.append(float(row[3]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
mu1=mu[0:575]
n=0.0001
pO1=np.multiply(pO,n)
k=0  
for i in range(6):
    for j in range(575):
        data1[5-i][j]=O[j+k]
        data4[5-i][j]=pO1[j+k]
        data5[5-i][j]=lm[j+k]
    k=k+575

ax1=plt.subplot(3, 3, 1)
ax1.text(-0.18, 1.05, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("density")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\rho$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,0.0125)
plt.plot(mu1, data4[0][:], c = 'brown')
plt.plot(mu1, data4[1][:], c = 'r')
plt.plot(mu1, data4[2][:], c = 'g')
plt.plot(mu1, data4[3][:], c = 'c')
plt.plot(mu1, data4[4][:], c = 'b')
plt.plot(mu1, data4[5][:], c = 'k')
plt.axhline(y = 0.011, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 0.0068, color = 'grey', linestyle = '-') 
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")

mux=[]
E0x=[]
Ox=[]
pOx=[]
lmx=[]
data1x = np.zeros( (6, 575) )
data4x = np.zeros( (6, 575) )
data5x = np.zeros( (6, 575) )
f = open('/lnZ_nondegen/fixed_tonks_gas_model/density_vs_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mux.append(float(row[0]))
    E0x.append(float(row[1]))
    Ox.append(float(row[3]))
    pOx.append(float(row[6]))
    lmx.append(float(row[7]))
f.close()
mu1x=mux[0:575]
n=0.0001
pO1x=np.multiply(pOx,n)
k=0  
for i in range(6):
    for j in range(575):
        data1x[5-i][j]=Ox[j+k]
        data4x[5-i][j]=pO1x[j+k]
        data5x[5-i][j]=lmx[j+k]
    k=k+575
         
plt.plot(mu1x, data4x[0][:], c = 'brown',linestyle='dashed')
plt.plot(mu1x, data4x[1][:], c = 'r',linestyle='dashed')
plt.plot(mu1x, data4x[2][:], c = 'g',linestyle='dashed')
plt.plot(mu1x, data4x[3][:], c = 'c',linestyle='dashed')
plt.plot(mu1x, data4x[4][:], c = 'b',linestyle='dashed')
plt.plot(mu1x, data4x[5][:], c = 'k',linestyle='dashed')
print('Panel 1.......')
#------------------------------------------------------------------

ax1=plt.subplot(3, 3, 2)
ax1.text(-0.18, 1.05, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Mean Occupancy")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('$O$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,1)
plt.plot(mu1, data1[0][:], c = 'brown')
plt.plot(mu1, data1[1][:], c = 'r')
plt.plot(mu1, data1[2][:], c = 'g')
plt.plot(mu1, data1[3][:], c = 'c')
plt.plot(mu1, data1[4][:], c = 'b')
plt.plot(mu1, data1[5][:], c = 'k')
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")

plt.plot(mu1x, data1x[0][:], c = 'brown',linestyle='dashed')
plt.plot(mu1x, data1x[1][:], c = 'r',linestyle='dashed')
plt.plot(mu1x, data1x[2][:], c = 'g',linestyle='dashed')
plt.plot(mu1x, data1x[3][:], c = 'c',linestyle='dashed')
plt.plot(mu1x, data1x[4][:], c = 'b',linestyle='dashed')
plt.plot(mu1x, data1x[5][:], c = 'k',linestyle='dashed')
print('Panel 2.......')

ax1=plt.subplot(3, 3, 3)
ax1.text(-0.18, 1.05, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Mean length")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('$\\langle l\\rangle$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(50,160)
plt.plot(mu1, data5[0][:], c = 'brown')
plt.plot(mu1, data5[1][:], c = 'r')
plt.plot(mu1, data5[2][:], c = 'g')
plt.plot(mu1, data5[3][:], c = 'c')
plt.plot(mu1, data5[4][:], c = 'b')
plt.plot(mu1, data5[5][:], c = 'k')
plt.axhline(y = 91, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 147, color = 'grey', linestyle = '-') 
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 3.......')

mu=[]
E0=[]
O=[]
pO=[]
lm=[]
data1 = np.zeros( (6, 575) )
data4 = np.zeros( (6, 575) )
data5 = np.zeros( (6, 575) )
f = open('/lnZ_nondegen/no_nieb_tomchou/density_amin91_vs_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    O.append(float(row[3]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
mu1=mu[0:575]
n=0.0001
pO1=np.multiply(pO,n)
k=0     
for i in range(6):
    for j in range(575):
        data1[5-i][j]=O[j+k]
        data4[5-i][j]=pO1[j+k]
        data5[5-i][j]=lm[j+k]
    k=k+575
         
#--------theory+sim-----------------------------------
mu=[]
lm=[]
data3 = np.zeros( (6, 575) )
f = open('/lnZ_nondegen/fixed_tonks_gas_model/density_vs_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    lm.append(float(row[7]))
f.close()
mu1=mu[0:575]
n=0.0001*20
pO1=np.multiply(pO,n)
k=0  
for i in range(6):
    for j in range(575):
        data3[5-i][j]=lm[j+k]
    k=k+575
    
#axins = ax1.inset_axes([0.60, 0.65, 0.5, 0.52]) #[left, bottom, width, height]

plt.plot(mu1, data3[0][:], c = 'brown',linestyle='dashed')
plt.plot(mu1, data3[1][:], c = 'r',linestyle='dashed')
plt.plot(mu1, data3[2][:], c = 'g',linestyle='dashed')
plt.plot(mu1, data3[3][:], c = 'c',linestyle='dashed')
plt.plot(mu1, data3[4][:], c = 'b',linestyle='dashed')
plt.plot(mu1, data3[5][:], c = 'k',linestyle='dashed')
#plt.set_xlim(-75, 40)
#plt.set_ylim(50,160)


ax1=plt.subplot(3, 3, 4)
ax1.text(-0.18, 1.05, string.ascii_uppercase[3], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("density")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\rho$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,0.0125)
plt.plot(mu1, data4[0][:], c = 'brown')
plt.plot(mu1, data4[1][:], c = 'r')
plt.plot(mu1, data4[2][:], c = 'g')
plt.plot(mu1, data4[3][:], c = 'c')
plt.plot(mu1, data4[4][:], c = 'b')
plt.plot(mu1, data4[5][:], c = 'k')
plt.axhline(y = 0.009, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 0.006, color = 'grey', linestyle = '-') 
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 4.......')

ax1=plt.subplot(3, 3, 5)
ax1.text(-0.18, 1.05, string.ascii_uppercase[4], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Mean Occupancy")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('$O$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,1)
plt.plot(mu1, data1[0][:], c = 'brown')
plt.plot(mu1, data1[1][:], c = 'r')
plt.plot(mu1, data1[2][:], c = 'g')
plt.plot(mu1, data1[3][:], c = 'c')
plt.plot(mu1, data1[4][:], c = 'b')
plt.plot(mu1, data1[5][:], c = 'k')
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 5.......')

ax1=plt.subplot(3, 3, 6)
ax1.text(-0.18, 1.05, string.ascii_uppercase[5], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Mean Occupancy")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('$\\langle l\\rangle$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(50,160)
plt.plot(mu1, data5[0][:], c = 'brown')
plt.plot(mu1, data5[1][:], c = 'r')
plt.plot(mu1, data5[2][:], c = 'g')
plt.plot(mu1, data5[3][:], c = 'c')
plt.plot(mu1, data5[4][:], c = 'b')
plt.plot(mu1, data5[5][:], c = 'k')
plt.axhline(y = 91, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 147, color = 'grey', linestyle = '-') 
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 6.......')

mu=[]
E0=[]
O=[]
pO=[]
lm=[]
data1 = np.zeros( (6, 575) )
data4 = np.zeros( (6, 575) )
data5 = np.zeros( (6, 575) )
f = open('/lnZ_nondegen/no_nieb_tomchou/density_amin1_vs_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    O.append(float(row[3]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
mu1=mu[0:575]
n=0.0001
pO1=np.multiply(pO,n)
k=0   
for i in range(6):
    for j in range(575):
        data1[5-i][j]=O[j+k]
        data4[5-i][j]=pO1[j+k]
        data5[5-i][j]=lm[j+k]
    k=k+575         

ax1=plt.subplot(3, 3, 7)
ax1.text(-0.18, 1.05, string.ascii_uppercase[6], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("density")
plt.xlabel(r'$\mu$',fontsize=12)
plt.ylabel(r'$\rho$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,0.05)
plt.plot(mu1, data4[0][:], c = 'brown',label=r'$\epsilon_{0}=-0.411$')
plt.plot(mu1, data4[1][:], c = 'r',label=r'$\epsilon_{0}=-0.274$')
plt.plot(mu1, data4[2][:], c = 'g',label=r'$\epsilon_{0}=-0.137$')
plt.plot(mu1, data4[3][:], c = 'c',label=r'$\epsilon_{0}=-0.034$')
plt.plot(mu1, data4[4][:], c = 'b',label=r'$\epsilon_{0}=0$')
plt.plot(mu1, data4[5][:], c = 'k',label=r'$\epsilon_{0}=0.068$')
plt.legend(loc="upper left",bbox_to_anchor=(-0.02, 0.9),borderpad=0.25,labelspacing=0.25,fontsize=10,frameon=False)
plt.axhline(y = 0.048, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 0.006, color = 'grey', linestyle = '-') 
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 7.......')

ax1=plt.subplot(3, 3, 8)
ax1.text(-0.18, 1.05, string.ascii_uppercase[7], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Mean Occupancy")
plt.xlabel(r'$\mu$',fontsize=12)
plt.ylabel('$O$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,1)
plt.plot(mu1, data1[0][:], c = 'brown')
plt.plot(mu1, data1[1][:], c = 'r')
plt.plot(mu1, data1[2][:], c = 'g')
plt.plot(mu1, data1[3][:], c = 'c')
plt.plot(mu1, data1[4][:], c = 'b')
plt.plot(mu1, data1[5][:], c = 'k')
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 8.......')

ax1=plt.subplot(3, 3, 9)
ax1.text(-0.18, 1.05, string.ascii_uppercase[8], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Mean length")
plt.xlabel(r'$\mu$',fontsize=12)
plt.ylabel('$\\langle l\\rangle$',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,200)
plt.plot(mu1, data5[0][:], c = 'brown')
plt.plot(mu1, data5[1][:], c = 'r')
plt.plot(mu1, data5[2][:], c = 'g')
plt.plot(mu1, data5[3][:], c = 'c')
plt.plot(mu1, data5[4][:], c = 'b')
plt.plot(mu1, data5[5][:], c = 'k')
plt.axhline(y = 1, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 147, color = 'grey', linestyle = '-') 
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
print('Panel 9.......')
    
fig.tight_layout()
plt.subplots_adjust(wspace=0.29, hspace=0.25)
plt.show()
  
 




