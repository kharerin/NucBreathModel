import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import string


mu1=[]
E01=[]
f = open('occ_0p825.txt','r')
for row in f:
    row = row.split(' ')
    mu1.append(float(row[0]))
    E01.append(float(row[1]))
f.close()  
n=-56
E11=np.multiply(E01,n)
#mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,Osum3/L1,Osum2/L1,Osum1/L1  

mu=[]
E0=[]
dO=[]
O=[]
rms=[]
rms2=[]
pO=[]
lm=[]
data1 = np.zeros( (500, 500) )
data2 = np.zeros( (500, 500) )
data3 = np.zeros( (500, 500) )
data4 = np.zeros( (500, 500) )
data5 = np.zeros( (500, 500) )
f = open('d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    dO.append(float(row[2]))
    O.append(float(row[3]))
    rms.append(float(row[4]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
n=-56
E1=np.multiply(E0,n) 
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)
k=0   
for i in range(500):
    for j in range(500):
        data1[499-i][j]=O[j+k]
        data2[499-i][j]=dO[j+k]
        data3[499-i][j]=rms21[j+k]
        data4[499-i][j]=pO1[j+k]
        data5[499-i][j]=lm[j+k]
    k=k+500
            
fig, ax = plt.subplots(2,4,figsize=(15,8))
#fig, ax = plt.subplots(constrained_layout=True)

ax1=plt.subplot(2, 4, 1)
ax1.text(-0.1, 1.1, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')

plt.title("Mean occupancy O")
plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('$E_{0d}$',fontsize=12)

plt.xlim(-60, 40)
plt.ylim(25,-25)

imx1=plt.imshow(data1,extent=[-60, 40,25,-25], vmin=0, vmax=1, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.plot(mu1,E11, linewidth=1, c = 'k')

ax1=plt.subplot(2, 4, 2)
ax1.text(-0.1, 1.1, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
plt.title("$\Delta$O",fontsize=12)
plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
imx2=plt.imshow(data2,extent=[-60, 40,25,-25], vmin=data2.min(), vmax=data2.max(), cmap='jet',aspect=2);
plt.colorbar(imx2,shrink=0.5)

ax1=plt.subplot(2, 4, 3)
ax1.text(-0.1, 1.1, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r'Density ($\rho$)',fontsize=12)
plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
imx4=plt.imshow(data4,extent=[-60, 40,25,-25], vmin=0, vmax=0.012, cmap='jet',aspect=2);
plt.colorbar(imx4,shrink=0.5)

ax1=plt.subplot(2, 4, 4)
ax1.text(-0.1, 1.1, string.ascii_uppercase[3], transform=ax1.transAxes, size=20, weight='bold')
plt.title("Rmsd",fontsize=12)
plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
imx3=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=data3.min(), vmax=data3.max(), cmap='jet',aspect=2);
plt.colorbar(imx3,shrink=0.5)

ax1=plt.subplot(2, 4, 5)
ax1.text(-0.1, 1.1, string.ascii_uppercase[3], transform=ax1.transAxes, size=20, weight='bold')
plt.title("mean_l",fontsize=12)
plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
imx3=plt.imshow(data5,extent=[-60, 40,25,-25], vmin=1, vmax=150, cmap='jet',aspect=2);
plt.colorbar(imx3,shrink=0.5)

mu=[]
E0=[]
dO=[]
rms=[]
rms2=[]
pO=[]
lm=[]
O1=[]
O2=[]
O3=[]

f = open('occ_0p825.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    dO.append(float(row[2]))
    rms.append(float(row[4]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
    O3.append(float(row[8]))
    O2.append(float(row[9]))
    O1.append(float(row[10]))
f.close()
n=-56
E1=np.multiply(E0,n) 
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)  

ax1=plt.subplot(2, 4, 6)
ax1.text(-0.1, 1.1, string.ascii_uppercase[4], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("$\Delta$O")
plt.xlabel('$E_{0d}$',fontsize=12)
plt.ylabel('$\Delta$O',fontsize=12)
plt.xlim(-25, 25)
plt.ylim(-0.4,0.4)
plt.plot(E1, dO, c = 'k')
plt.vlines(-2.65, -0.4, 0.4, linestyles ="dotted", colors ="k")

ax1=plt.subplot(2, 4, 7)
ax1.text(-0.1, 1.1, string.ascii_uppercase[5], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("Density")
plt.xlabel('$E_{0d}$',fontsize=12)
plt.ylabel(r'$\rho$',fontsize=12)
plt.xlim(-25, 25)
plt.ylim(0,0.012)
plt.plot(E1, pO1, c = 'k')
plt.vlines(-2.65, 0, 0.012, linestyles ="dotted", colors ="k")


    
#fig.tight_layout()
plt.show()
  
 




