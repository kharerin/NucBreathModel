import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import math
import numpy as np
import string

fig, ax = plt.subplots(4,3,figsize=(9,10))

lmin=91
lmax=147
dle=lmax-lmin
mu1xx=[]
E01xx=[]
for j in range(500):
    h=-250+j
    E0=h/(10*(147-91)); 
    if(E0>0 or E0<0):
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    else:
       E0=0.00001;
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    mux=(-1)*E0*le-0.34-math.log(le) #mux=(-1)*E0*le+1-math.log(le)
    mu1xx.append(mux)
    E01xx.append(E0)
n=-56
E11xx=np.multiply(E01xx,n) 
#mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1  
mu2=[]
E02=[]
O=[]
rms2=[]
pO=[]
lm=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/lnZ_nondegen/mu_vs_E0d/mean_l/d1_mu_E0d.txt','r') # model 1
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    O.append(float(row[3]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close() 
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)

mu1=[]
E01=[]
mu1x=[]
E01x=[]
mu1xr=[]
E01xr=[]
mu1xO1=[]
E01xO1=[]
mu1xO2=[]
E01xO2=[]
mu1xO3=[]
E01xO3=[]
for i in range(250*250):
    if lm[i]>118.4 and lm[i]<119.6:
       mu1x.append(mu2[i])
       E01x.append(E02[i])
    if rms21[i]>0.000001 and rms21[i]<0.000005 and O[i]>0.4:
       mu1xr.append(mu2[i])
       E01xr.append(E02[i])
    if O[i]>0.799 and O[i]<0.801:
       mu1.append(mu2[i])
       E01.append(E02[i])
    if O[i]>0.849 and O[i]<0.851:
       mu1xO1.append(mu2[i])
       E01xO1.append(E02[i])
    if O[i]>0.899 and O[i]<0.901:
       mu1xO2.append(mu2[i])
       E01xO2.append(E02[i])
    if O[i]>0.497 and O[i]<0.503:
       mu1xO3.append(mu2[i])
       E01xO3.append(E02[i])
    if (lm[i]>120 and lm[i]<130) and (O[i]>0.7 and O[i]<0.9) and (pO1[i]>0.006 and pO1[i]<0.007): 
       #mu1xr.append(mu2[i])
       #E01xr.append(E02[i])
       print('>>>>>>>>>1st model...i...mu...e0...occ',mu2[i],E02[i],O[i])
       
n=-56
E11x=np.multiply(E01x,n)  
E11xr=np.multiply(E01xr,n)   
E11=np.multiply(E01,n)
E11xO1=np.multiply(E01xO1,n)
E11xO2=np.multiply(E01xO2,n) 
E11xO3=np.multiply(E01xO3,n) 
mu1y=[]
E01y=[]
pO1y=[]
f = open('/lnZ_nondegen/no_nieb_cedric/d1_mu_E0d.txt','r') # model 1 & theory 
for row in f:
    row = row.split(' ')
    mu1y.append(float(row[0]))
    E01y.append(float(row[1]))
    pO1y.append(float(row[2]))
f.close()  
n=56
E11y=np.multiply(E01y,n)
k=0   
for i in range(250):
    for j in range(250):
        data1[249-i][j]=O[j+k]
        data2[249-i][j]=pO1[j+k]
        data3[249-i][j]=lm[j+k]
    k=k+250
            
ax1=plt.subplot(4, 3, 1)
ax1.text(-0.16, 1.0, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
ax1.text(0.3, 0.8, 'ID', transform=ax1.transAxes, size=14, weight='bold', color='black')
ax1.text(0.75, 0.6, 'HD', transform=ax1.transAxes, size=14, weight='bold', color='black')
ax1.text(0.2, 0.3, 'LD', transform=ax1.transAxes, size=14, weight='bold', color='white')
plt.title(r'Density $\rho$',fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_{0}$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data2,extent=[-60, 40,25,-25], vmin=0, vmax=0.012, cmap='jet',aspect=2);
#plt.setp(ax1.get_xticklabels(), visible=False)
#plt.setp(ax1.get_yticklabels(), visible=False)
plt.xticks([])
#plt.yticks([])
plt.colorbar(imx1,shrink=0.5)
plt.plot(mu1y,E11y, linewidth=1, c = 'k')
plt.plot(mu1x,E11x, linestyle='dotted', linewidth=1.5, c = 'w')
plt.plot(mu1xx,E11xx, linewidth=1.5, c = 'r',zorder=0, label=r'$\mu$*') #linestyle='dashed'
plt.legend(loc="lower left",fontsize=9,frameon=True,facecolor='white') #bbox_to_anchor=(1.4, -0.1)
plt.scatter(-12.8,-0.057*56, marker = '+', s=50, c = 'k',zorder=1)

ax1=plt.subplot(4, 3, 2)
ax1.text(-0.16, 1.0, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
plt.title("Occupancy $O$")
#plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('$E_{0d}$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)

imx1=plt.imshow(data1,extent=[-60, 40,25,-25], vmin=0, vmax=1, cmap='jet',aspect=2);
#plt.setp(ax1.get_xticklabels(), visible=False)
plt.setp(ax1.get_yticklabels(), visible=False)
plt.xticks([])
plt.yticks([])
plt.colorbar(imx1,shrink=0.5)
plt.plot(mu1xO2,E11xO2, linewidth=1, linestyle='dashed', c = 'k',label='90%')
plt.plot(mu1xO1,E11xO1, linewidth=1, linestyle='dotted', c = 'k',label='85%')
plt.plot(mu1,E11, linewidth=1, c = 'k',label='80%')
plt.plot(mu1xO3,E11xO3, linewidth=1, c = 'r',linestyle='dashed',label='50%')
plt.legend(loc="lower left",fontsize=9,frameon=True) #bbox_to_anchor=(1.43, -0.28)
plt.scatter(-12.8,-0.057*56, marker = '+', s=50, c = 'k',zorder=1)

ax1=plt.subplot(4, 3, 3)
ax1.text(-0.16, 1.0, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"Wrapping length $\langle l \rangle$",fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)
imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=0, vmax=150, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.plot(mu1x,E11x, linewidth=1, c = 'k')
plt.xticks([])
plt.yticks([])
plt.scatter(-12.8,-0.057*56, marker = '+', s=50, c = 'k',zorder=1)

lmin=91
lmax=147
dle=lmax-lmin
mu1xx=[]
E01xx=[]
for j in range(500):
    h=-250+j
    E0=h/(10*(147-91)); 
    if(E0>0 or E0<0):
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    else:
       E0=0.00001;
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    mux=(-1)*E0*le-0.34-math.log(le+20)
    mu1xx.append(mux)
    E01xx.append(E0)
n=-56
E11xx=np.multiply(E01xx,n) 
#mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1 
mu2=[]
E02=[]
O=[]
rms2=[]
pO=[]
lm=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d/mean_l/d1_mu_E0d.txt','r') # model 2
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    O.append(float(row[3]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)

mu1=[]
E01=[]
mu1x=[]
E01x=[]
mu1xr=[]
E01xr=[]
mu1xO1=[]
E01xO1=[]
mu1xO2=[]
E01xO2=[]
mu1xO3=[]
E01xO3=[]
for i in range(250*250):
    if lm[i]>118.4 and lm[i]<119.6:
       mu1x.append(mu2[i])
       E01x.append(E02[i])
    if rms21[i]>0.000001 and rms21[i]<0.000005 and O[i]>0.4:
       mu1xr.append(mu2[i])
       E01xr.append(E02[i])
    if O[i]>0.799 and O[i]<0.801:
       mu1.append(mu2[i])
       E01.append(E02[i])
    if O[i]>0.699 and O[i]<0.701:
       mu1xO1.append(mu2[i])
       E01xO1.append(E02[i])
    if O[i]>0.749 and O[i]<0.751:
       mu1xO2.append(mu2[i])
       E01xO2.append(E02[i])  
    if O[i]>0.497 and O[i]<0.503:
       mu1xO3.append(mu2[i])
       E01xO3.append(E02[i])
    if (lm[i]>120 and lm[i]<130) and (O[i]>0.7 and O[i]<0.9) and (pO1[i]>0.006 and pO1[i]<0.007): 
       #mu1xr.append(mu2[i])
       #E01xr.append(E02[i])
       print('>>>>>>>>>2nd model...i...mu...e0...occ',mu2[i],E02[i],O[i])

n=-56
E11x=np.multiply(E01x,n) 
E11xr=np.multiply(E01xr,n) 
E11=np.multiply(E01,n)
E11xO1=np.multiply(E01xO1,n)
E11xO2=np.multiply(E01xO2,n) 
E11xO3=np.multiply(E01xO3,n)
k=0   
for i in range(250):
    for j in range(250):
        data1[249-i][j]=O[j+k]
        data2[249-i][j]=pO1[j+k]
        data3[249-i][j]=lm[j+k]
    k=k+250

ax1=plt.subplot(4, 3, 4)
ax1.text(-0.16, 1.0, string.ascii_uppercase[3], transform=ax1.transAxes, size=20, weight='bold')
ax1.text(0.3, 0.8, 'ID', transform=ax1.transAxes, size=14, weight='bold', color='black')
ax1.text(0.75, 0.6, 'HD', transform=ax1.transAxes, size=14, weight='bold', color='black')
ax1.text(0.2, 0.3, 'LD', transform=ax1.transAxes, size=14, weight='bold', color='white')
#plt.title(r'Density ($\rho$)',fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_{0}$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data2,extent=[-60, 40,25,-25], vmin=0, vmax=0.012, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.xticks([])
#plt.yticks([])
plt.plot(mu1xr,E11xr, linewidth=1, c = 'k')
plt.plot(mu1x,E11x, linestyle='dotted', linewidth=1.5, c = 'w')
plt.plot(mu1xx,E11xx, linewidth=1.5, c = 'r', zorder=0,label=r'$\mu$*') #linestyle='dashed'
plt.legend(loc="lower left",fontsize=9,frameon=True) #bbox_to_anchor=(1.4, -0.1)
plt.scatter(-9.2,-0.064*56, marker = '+', s=50, c = 'k',zorder=1)

ax1=plt.subplot(4, 3, 5)
ax1.text(-0.16, 1.0, string.ascii_uppercase[4], transform=ax1.transAxes, size=20, weight='bold')

#plt.title("Mean occupancy O")
#plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('$E_{0d}$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)

imx1=plt.imshow(data1,extent=[-60, 40,25,-25], vmin=0, vmax=1, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.xticks([])
plt.yticks([])

plt.plot(mu1,E11, linewidth=1, c = 'k',label='80%')
plt.plot(mu1xO2,E11xO2, linewidth=1, linestyle='dashed', c = 'k',label='75%')
plt.plot(mu1xO1,E11xO1, linewidth=1, linestyle='dotted', c = 'k',label='70%')
plt.plot(mu1xO3,E11xO3, linewidth=1, linestyle='dashed', c = 'r',label='50%')
plt.legend(loc="lower left",fontsize=9,frameon=True) #bbox_to_anchor=(1.53, -0.25)
plt.scatter(-9.2,-0.064*56, marker = '+', s=50, c = 'k',zorder=1)

ax1=plt.subplot(4, 3, 6)
ax1.text(-0.16, 1.0, string.ascii_uppercase[5], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("<l>",fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)
imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=0, vmax=150, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.plot(mu1x,E11x, linewidth=1, c = 'k')
plt.xticks([])
plt.yticks([])
plt.scatter(-9.2,-0.064*56, marker = '+', s=50, c = 'k',zorder=1)

lmin=1
lmax=147
dle=lmax-lmin
mu1xx=[]
E01xx=[]
for j in range(500):
    h=-250+j
    E0=h/(10*(147-91)); 
    if(E0>0 or E0<0):
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    else:
       E0=0.00001;
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    mux=(-1)*E0*le-0.34-math.log(le+20)
    mu1xx.append(mux)
    E01xx.append(E0)
n=-56
E11xx=np.multiply(E01xx,n) 
#mu,E0,(Omd[56]-Omd[0]),Osum/L1,d1,d2,pO,lm,Osum3/L1,Osum2/L1,Osum1/L1  
mu2=[]
E02=[]
O=[]
rms2=[]
pO=[]
lm=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d_amin1/mean_l/d1_mu_E0d.txt','r') # model 3
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    O.append(float(row[3]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)

mu1=[]
E01=[]
mu1x=[]
E01x=[]
mu1xr=[]
E01xr=[]
mu1xO1=[]
E01xO1=[]
mu1xO2=[]
E01xO2=[]
mu1xO3=[]
E01xO3=[]
for i in range(250*250):
    if lm[i]>73 and lm[i]<75:
       mu1x.append(mu2[i])
       E01x.append(E02[i])
    if rms21[i]>0.000001 and rms21[i]<0.000005 and O[i]>0.4:
       mu1xr.append(mu2[i])
       E01xr.append(E02[i])
    if O[i]>0.699 and O[i]<0.701:
       mu1xO1.append(mu2[i])
       E01xO1.append(E02[i])
    if O[i]>0.749 and O[i]<0.751:
       mu1xO2.append(mu2[i])
       E01xO2.append(E02[i])
    if O[i]>0.497 and O[i]<0.503:
       mu1xO3.append(mu2[i])
       E01xO3.append(E02[i])
    if (lm[i]>120 and lm[i]<130) and (O[i]>0.7 and O[i]<0.9) and (pO1[i]>0.006 and pO1[i]<0.007): 
       #mu1xr.append(mu2[i])
       #E01xr.append(E02[i])
       print('>>>>>>>>>3rd model...i...mu...e0...occ',mu2[i],E02[i],O[i])
n=-56
E11x=np.multiply(E01x,n) 
E11xr=np.multiply(E01xr,n)
E11xO1=np.multiply(E01xO1,n)
E11xO2=np.multiply(E01xO2,n) 
E11xO3=np.multiply(E01xO3,n) 
k=0   
for i in range(250):
    for j in range(250):
        data1[249-i][j]=O[j+k]
        data2[249-i][j]=pO1[j+k]
        data3[249-i][j]=lm[j+k]
    k=k+250

ax1=plt.subplot(4, 3, 7)
ax1.text(-0.16, 1.0, string.ascii_uppercase[6], transform=ax1.transAxes, size=20, weight='bold')
ax1.text(0.3, 0.8, 'ID', transform=ax1.transAxes, size=14, weight='bold', color='black')
ax1.text(0.75, 0.6, 'HD', transform=ax1.transAxes, size=14, weight='bold', color='black')
ax1.text(0.2, 0.3, 'LD', transform=ax1.transAxes, size=14, weight='bold', color='white')
#plt.title(r'Density ($\rho$)',fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_{0}$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data2,extent=[-60, 40,25,-25], vmin=0, vmax=0.012, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.xticks([])
#plt.yticks([])
plt.plot(mu1xr,E11xr, linewidth=1, c = 'k')
plt.plot(mu1x,E11x, linestyle='dotted', linewidth=1.5, c = 'w')
plt.plot(mu1xx,E11xx, linewidth=1.5, c = 'r',zorder=0, label=r'$\mu$*') #linestyle='dashed'
plt.legend(loc="lower left",fontsize=9,frameon=True) #,bbox_to_anchor=(1.4, -0.1)
plt.scatter(-11.6,-0.096*56, marker = '+', s=50, c = 'k',zorder=1)

ax1=plt.subplot(4, 3, 8)
ax1.text(-0.16, 1.0, string.ascii_uppercase[7], transform=ax1.transAxes, size=20, weight='bold')

#plt.title("Mean occupancy O")
#plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('$E_{0d}$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)

imx1=plt.imshow(data1,extent=[-60, 40,25,-25], vmin=0, vmax=1, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.xticks([])
plt.yticks([])
#plt.plot(mu1,E11, linewidth=1, c = 'k',label='80%')
#plt.plot(mu1xO2,E11xO2, linewidth=1, linestyle='dashed', c = 'k',label='75%')
model3 = np.poly1d(np.polyfit(mu1xO3, E11xO3, 12))
polyline = np.linspace(-60, 10, 100)
#plt.scatter(mu1xO3,E11xO3, s=5, linestyle='dashed', c = 'r')
plt.plot(polyline, model3(polyline), c='r',linewidth=1,linestyle='dashed',label='50%')
plt.legend(loc="lower left",fontsize=9,frameon=True) #bbox_to_anchor=(1.53, -0.2)
plt.scatter(-11.6,-0.096*56, marker = '+', s=50, c = 'k',zorder=1)

ax1=plt.subplot(4, 3, 9)
ax1.text(-0.16, 1.0, string.ascii_uppercase[8], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("<l>",fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
#plt.ylabel('E0d')
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)
imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=0, vmax=150, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
plt.plot(mu1x,E11x, linewidth=1, c = 'k')
plt.xticks([])
plt.yticks([])
plt.scatter(-11.6,-0.096*56, marker = '+', s=50, c = 'k',zorder=1)

#--------theory+sim--------------------------------------
print('>>>>>>>>>theory+sim...end!')
lmin=91
lmax=147
dle=lmax-lmin
mu1xx=[]
E01xx=[]
for j in range(500):
    h=-250+j
    E0=h/(10*(147-91)); 
    if(E0>0 or E0<0):
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    else:
       E0=0.00001;
       le=lmin+dle*(math.exp(-(-1)*E0*dle)/(math.exp(-(-1)*E0*dle)-1))+1/((-1)*E0)
    mux=(-1)*E0*le-0.34-math.log(le)
    mu1xx.append(mux)
    E01xx.append(E0)
n=-56
E11xx=np.multiply(E01xx,n) 
 
mu2=[]
E02=[]
O=[]
rms2=[]
pO=[]
lm=[]
data1 = np.zeros( (500, 500) )
data2 = np.zeros( (500, 500) )
data3 = np.zeros( (500, 500) )
f = open('/lnZ_nondegen/fixed_tonks_gas_model/mu_vs_E0d/d1_mu_E0d.txt','r') # fixed-size model
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    O.append(float(row[3]))
    rms2.append(float(row[5]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close() 
n=0.0001
rms21=np.multiply(rms2,n)
pO1=np.multiply(pO,n)

mu1x=[]
E01x=[]
mu1xr=[]
E01xr=[]
mu1xO1=[]
E01xO1=[]
mu1xO2=[]
E01xO2=[]
for i in range(500*500):
    if O[i]>0.799 and O[i]<0.801:
       mu1xO1.append(mu2[i])
       E01xO1.append(E02[i])
    if O[i]>0.498 and O[i]<0.502:
       mu1xO2.append(mu2[i])
       E01xO2.append(E02[i])
       #print(mu2[i],E02[i],O[i])

n=-56
E11x=np.multiply(E01x,n)  
E11xr=np.multiply(E01xr,n)   
E11xO1=np.multiply(E01xO1,n)
E11xO2=np.multiply(E01xO2,n) 
k=0   
for i in range(500):
    for j in range(500):
        data1[499-i][j]=pO1[j+k]
        data2[499-i][j]=O[j+k]
        data3[499-i][j]=lm[j+k]
    k=k+500  

ax1=plt.subplot(4, 3, 10)
ax1.text(-0.16, 1.0, string.ascii_uppercase[9], transform=ax1.transAxes, size=20, weight='bold')
plt.ylabel(r'$\epsilon_{0}$',fontsize=12)
plt.xlabel(r'$\mu$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data1,extent=[-60, 40,25,-25], vmin=0, vmax=0.012, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)  
#plt.xticks([])
#plt.yticks([])
plt.plot(mu1xx,E11xx, linewidth=1.5, c = 'r', label=r'$\mu$*') #linestyle='dashed'
plt.legend(loc="lower left",fontsize=9,frameon=True) #bbox_to_anchor=(1.4, -0.1)

ax1=plt.subplot(4, 3, 11)
ax1.text(-0.16, 1.0, string.ascii_uppercase[10], transform=ax1.transAxes, size=20, weight='bold')
#plt.ylabel(r'Density ($\rho$)',fontsize=12)
plt.xlabel(r'$\mu$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)

imx1=plt.imshow(data2,extent=[-60, 40,25,-25], vmin=0, vmax=1, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5) 
#plt.xticks([])
plt.yticks([])
plt.plot(mu1xO1,E11xO1, linewidth=1, c = 'k',label='80%') #linestyle='dashed'
plt.plot(mu1xO2,E11xO2, linewidth=1, c = 'r',linestyle='dashed',label='50%') #linestyle='dashed'
plt.legend(loc="lower left",fontsize=9,frameon=True) #bbox_to_anchor=(1.53, -0.2)

ax1=plt.subplot(4, 3, 12)
ax1.text(-0.16, 1.0, string.ascii_uppercase[11], transform=ax1.transAxes, size=20, weight='bold')
#plt.ylabel(r'Density ($\rho$)',fontsize=12)
plt.xlabel(r'$\mu$',fontsize=12)
plt.xticks(np.arange(-50, 40, 25))
#plt.yticks(np.arange(90, 147, 10)) 
plt.xlim(-60, 40)
plt.ylim(25,-25)

imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=0, vmax=150, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5) 
#plt.xticks([])
plt.yticks([])

fig.tight_layout()
plt.subplots_adjust(wspace=0.15, hspace=0.08)
plt.show()
  
 




