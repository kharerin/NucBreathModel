import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import string

fig, ax = plt.subplots(1,2,figsize=(8.2,3.8))

mu2=[]
E02=[]
lm=[]
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_vs_E0d/mean_l/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    lm.append(float(row[7]))
f.close() 
mu1x=[]
E01x=[]
for i in range(250*250):
    if lm[i]>118.4 and lm[i]<119.6:
       mu1x.append(mu2[i])
       E01x.append(E02[i])
n=-56
E11x=np.multiply(E01x,n) 
mu=[]
E0=[]
dO=[]
O=[]
rms=[]
rms2=[]
pO=[]
data1 = np.zeros( (500, 500) )
data2 = np.zeros( (500, 500) )
data3 = np.zeros( (500, 500) )
data4 = np.zeros( (500, 500) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_vs_E0d/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    dO.append(float(row[2]))
    O.append(float(row[3]))
    rms.append(float(row[4]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
f.close()
n=-56
E1=np.multiply(E0,n) 
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)
mu1xr=[]
E01xr=[]
pO1xr=[]
for i in range(500*500):
    if rms21[i]>0.000001 and rms21[i]<0.000005 and O[i]>0.4:
       mu1xr.append(mu[i])
       E01xr.append(E0[i])
       pO1xr.append(pO1[i])

n=-56 
E11xr=np.multiply(E01xr,n)  
k=0   
for i in range(500):
    for j in range(500):
        data1[499-i][j]=O[j+k]
        data2[499-i][j]=dO[j+k]
        data3[499-i][j]=rms21[j+k]
        data4[499-i][j]=pO1[j+k]
    k=k+500

mu1=[]
E01=[]
pO1=[]
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_cedric/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu1.append(float(row[0]))
    E01.append(float(row[1]))
    pO1.append(float(row[2]))
f.close()  
n=56
E11=np.multiply(E01,n)
n=-56
E12=np.multiply(E01,n)

ax1=plt.subplot(1, 2, 1)
ax1.text(-0.2, 1.1, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\rm{RMSD}_{\rho}$",fontsize=12)
plt.xlabel(r'$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_0$',fontsize=12)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 
imx3=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=data3.min(), vmax=data3.max(), cmap='jet',aspect=2);
plt.colorbar(imx3,shrink=0.5)
plt.plot(mu1xr,E11xr, linewidth=3, c = 'w',label=r'$\rm{RMSD}_{\rho}=0$') # zeros rmsd
plt.plot(mu1x,E11x, linewidth=1, c = 'magenta',label='$\\langle l \\rangle=119$') # iso <l>
plt.plot(mu1,E11, linewidth=1, c = 'k',linestyle='dashed',label='Analytical')  # analytic
plt.legend(loc="lower left",fontsize=10)

ax1=plt.subplot(1, 2, 2)
ax1.text(-0.2, 1.1, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
#plt.title('Critical line',fontsize=12)
plt.xlabel(r'$\epsilon_0$',fontsize=12)
plt.ylabel(r'$\rho$',fontsize=12)
E13=[0.001,0.0025,0.004,0.0055,0.007,0.0085,0.01]
labels = ['0.001','0.0025','0.004','0.0055','0.007','0.0085','0.01']
ax1.set_xticks([-25,-20, -15,-10, -5,0]) 
ax1.set_xticklabels(['$-0.45$', '$-0.36$', '$-0.27$', '$-0.18$', '$-0.09$', '0']) 
plt.yticks(E13, labels)
#plt.xlim(-60, 40)
plt.ylim(0.001,0.01)
#sc=plt.scatter(mu1xr, -E11xr, s=20, c=pO1xr,marker='o', vmin=0, vmax=0.012, cmap='jet')
#sc=plt.scatter(mu1, E12, s=20, c=pO1,marker='o',label='Analytic',vmin=0,vmax=0.012, cmap='jet') #facecolors='none', edgecolors='grey'
plt.plot(E11xr,pO1xr, linewidth=2, c = 'grey',linestyle='-',label='Model')  # analytic
plt.plot(E11,pO1, linewidth=1.5, c = 'k',linestyle='dashed',label='Analytical')  # analytic
#imx3=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=data3.min(), vmax=data3.max(), cmap='jet',aspect=2);
#plt.colorbar(sc,shrink=0.5)
print(min(pO1))
print(max(pO1))
print('Ok!....')
plt.legend(loc="lower left",fontsize=10)
    
fig.tight_layout()
fig.subplots_adjust(wspace=0.35, hspace=0.25, bottom=0.17, top=0.805)
plt.show()
  
 




