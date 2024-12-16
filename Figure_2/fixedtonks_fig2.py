import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np 
import math
import string

fig, ax = plt.subplots(figsize=(9,3.2))
#fig, ax = plt.subplots(constrained_layout=True)

plt.clf()

lmin=91
lmax=147
l = np.zeros(9)
p = np.zeros(10000)
o = np.zeros(10000)
ur = np.zeros(10000)
#cx = np.zeros(9)
l[0]=70; l[1]=80; l[2]=90; l[3]=100; l[4]=110; l[5]=120; l[6]=130; l[7]=140; l[8]=150; 
cx='rbgycmpkk'
print(cx[0])
ax=plt.subplot(1, 3, 1)
ax.text(-0.25, 1.07, string.ascii_uppercase[0], transform=ax.transAxes, size=20, weight='bold')
plt.ylabel(r'Density ($\rho$)',fontsize=12)
plt.xlabel(r'$\mu_r$',fontsize=12)
plt.xlim(-15, 15)
plt.ylim(0,0.015)
for j in range(9):
    px=[]
    ux=[]
    for i in range(10000):
        p[i]=i*0.00001+0.0000001
        o[i]=p[i]*l[j]
        #print(p[i],p[i]*l[j])
        if o[i]<1:
           ur[i]=math.log(p[i]*l[j])-math.log(1-p[i]*l[j])+p[i]*l[j]/(1-p[i]*l[j])-math.log(l[j])
           px.append(p[i])
           ux.append(ur[i])
    if len(ux)>10:
       plt.plot(ux, px) #linestyle='dotted'
       print('j.....',j,cx[j])
#plt.plot(ux, px, c='k') #linestyle='dotted' 

mu=[]
E0=[]
O=[]
pO=[]
lm=[]
data1 = np.zeros( (9, 575) )
data4 = np.zeros( (9, 575) )
data5 = np.zeros( (9, 575) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/density_fixedsize_vs_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    O.append(float(row[3]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
f.close()
x=np.arange(0, 575, 10)
mu1=mu[0:575]
n=0.0001
pO1=np.multiply(pO,n)
k=0  
for i in range(9):
    for j in range(575):
        data1[i][j]=O[j+k]
        data4[i][j]=pO1[j+k]
        data5[i][j]=lm[j+k]
    k=k+575

plt.plot(mu1, data4[0][:], linestyle='', c = '#1f77b4',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[1][:], linestyle='', c = '#ff7f0e',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[2][:], linestyle='', c = '#2ca02c',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[3][:], linestyle='', c = '#d62728',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[4][:], linestyle='', c = '#9467bd',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[5][:], linestyle='', c = '#8c564b',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[6][:], linestyle='', c = '#e377c2',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[7][:], linestyle='', c = '#7f7f7f',marker='o',ms=4,markevery=8)
plt.plot(mu1, data4[8][:], linestyle='', c = '#bcbd22',marker='o',ms=4,markevery=8)
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")
     

p = np.zeros(10000)
o = np.zeros(10000)
ur = np.zeros(10000)    
ax1=plt.subplot(1, 3, 2)
ax1.text(-0.25, 1.07, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
plt.ylabel('Occupancy ($O$)',fontsize=12)
plt.xlabel(r'$\mu_r$',fontsize=12)
plt.xlim(-15, 15)
plt.ylim(0,1)
for j in range(9):
    px=[]
    ux=[]
    ox=[]
    for i in range(10000):
        p[i]=i*0.00001+0.0000001
        o[i]=p[i]*l[j]
        #print(p[i],p[i]*l[j])
        if o[i]<1:
           ur[i]=math.log(p[i]*l[j])-math.log(1-p[i]*l[j])+p[i]*l[j]/(1-p[i]*l[j])-math.log(l[j])
           px.append(p[i])
           ux.append(ur[i])
           ox.append(o[i])
    if len(ux)>10:
       plt.plot(ux, ox,label='%dbp'%int(l[j])) #linestyle='dotted'
       print('j.....',j,cx[j])
#plt.plot(ux, ox, c='k') #linestyle='dotted' 

plt.plot(mu1, data1[0][:], linestyle='', c = '#1f77b4',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[1][:], linestyle='', c = '#ff7f0e',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[2][:], linestyle='', c = '#2ca02c',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[3][:], linestyle='', c = '#d62728',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[4][:], linestyle='', c = '#9467bd',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[5][:], linestyle='', c = '#8c564b',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[6][:], linestyle='', c = '#e377c2',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[7][:], linestyle='', c = '#7f7f7f',marker='o',ms=4,markevery=8)
plt.plot(mu1, data1[8][:], linestyle='', c = '#bcbd22',marker='o',ms=4,markevery=8)
#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")   
plt.legend(loc="upper right",bbox_to_anchor=(1.03, 0.85),borderpad=0.25,labelspacing=0.25,fontsize=10,frameon=False) 

mu=[]
E0=[]
O=[]
pO=[]
lm=[]
ms=[]
data1 = np.zeros(300)
data4 = np.zeros(300)
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/density_fixedsize_vs_mu/critical_mu/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    O.append(float(row[3]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
    ms.append(float(row[8]))
f.close() 
n=0.0001
pO1=np.multiply(pO,n)

muxx=[]
lsxx=[]
mux=np.zeros(7)
ls = np.zeros(7)
ls[0]=25; ls[1]=50; ls[2]=75; ls[3]=100; ls[4]=125; ls[5]=150; ls[6]=175; 
k=0  
for i in range(7):
    data1 = np.zeros(300)
    data5 = np.zeros(300)
    for j in range(300):
        data1[j]=mu[j+k]
        data5[j]=pO1[j+k]
        if pO1[j+k]>((0.31/ls[i])-0.0001) and pO1[j+k]<((0.31/ls[i])+0.0005) and ls[i]==25:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
        if pO1[j+k]>((0.31/ls[i])-0.0001) and pO1[j+k]<((0.31/ls[i])+0.0001) and ls[i]==50:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
        if pO1[j+k]>((0.31/ls[i])-0.0001) and pO1[j+k]<((0.31/ls[i])+0.0002) and ls[i]==75:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
        if pO1[j+k]>((0.31/ls[i])-0.0001) and pO1[j+k]<((0.31/ls[i])+0.0001) and ls[i]==100:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
        if pO1[j+k]>((0.31/ls[i])-0.00004) and pO1[j+k]<((0.31/ls[i])+0.00004) and ls[i]==125:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
        if pO1[j+k]>((0.31/ls[i])-0.00004) and pO1[j+k]<((0.31/ls[i])+0.00004) and ls[i]==150:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
        if pO1[j+k]>((0.31/ls[i])-0.00004) and pO1[j+k]<((0.31/ls[i])+0.00003) and ls[i]==175:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx.append(mu[j+k])
           lsxx.append(ls[i])
    #imx=list(data5).index(max(list(data5)))
    #mux[i]=data1[imx]
    #print('mx,imx..............',imx,data1[imx],max(list(data5)),ls[i])
    k=k+300


mu=[]
E0=[]
O=[]
pO=[]
lm=[]
ms=[]
data1 = np.zeros(300)
data4 = np.zeros(300)
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_-21p6kt_E0d8p4/nieb2/kevin_fig_ver2/large_remod2/score_scanA1/density_fixedsize_vs_mu/critical_mu_lk_20/den_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu.append(float(row[0]))
    E0.append(float(row[1]))
    O.append(float(row[3]))
    pO.append(float(row[6]))
    lm.append(float(row[7]))
    ms.append(float(row[8]))
f.close() 
n=0.0001
pO1=np.multiply(pO,n)

muxx1=[]
lsxx1=[]
mux=np.zeros(7)
ls = np.zeros(7)
ls[0]=25; ls[1]=50; ls[2]=75; ls[3]=100; ls[4]=125; ls[5]=150; ls[6]=175; 
k=0  
for i in range(7):
    data1 = np.zeros(300)
    data5 = np.zeros(300)
    for j in range(300):
        data1[j]=mu[j+k]
        data5[j]=pO1[j+k]
        if pO1[j+k]>((0.31/(ls[i]+20))-0.0001) and pO1[j+k]<((0.31/(ls[i]+20))+0.0005) and ls[i]==25:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
        if pO1[j+k]>((0.31/(ls[i]+20))-0.0001) and pO1[j+k]<((0.31/(ls[i]+20))+0.0001) and ls[i]==50:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
        if pO1[j+k]>((0.31/(ls[i]+20))-0.0001) and pO1[j+k]<((0.31/(ls[i]+20))+0.0001) and ls[i]==75:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
        if pO1[j+k]>((0.31/(ls[i]+20))-0.0001) and pO1[j+k]<((0.31/(ls[i]+20))+0.0001) and ls[i]==100:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
        if pO1[j+k]>((0.31/(ls[i]+20))-0.00004) and pO1[j+k]<((0.31/(ls[i]+20))+0.00004) and ls[i]==125:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
        if pO1[j+k]>((0.31/(ls[i]+20))-0.00004) and pO1[j+k]<((0.31/(ls[i]+20))+0.00004) and ls[i]==150:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
        if pO1[j+k]>((0.31/(ls[i]+20))-0.00004) and pO1[j+k]<((0.31/(ls[i]+20))+0.00003) and ls[i]==175:
           print('mx,imx..............',i,mu[j+k],ls[i])
           muxx1.append(mu[j+k])
           lsxx1.append(ls[i])
    #imx=list(data5).index(max(list(data5)))
    #mux[i]=data1[imx]
    #print('mx,imx..............',imx,data1[imx],max(list(data5)),ls[i])
    k=k+300

      
ax1=plt.subplot(1, 3, 3)
ax1.text(-0.25, 1.07, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
plt.ylabel(r'Chemical Potential ($\mu_r$*)',fontsize=12)
plt.xlabel('Wrapping length ($l$)',fontsize=12)
plt.xlim(0, 200)
plt.ylim(-8, 2)
l = np.zeros(20)
mu = np.zeros(20)
lk=0
for i in range(20):
    l[i]=i*10+5
    mu[i]=-0.34+math.log(1/(l[i]+lk))
    #print(p[i],p[i]*l[j])
plt.plot(l, mu,c='k',label='$l_{k}=0$') #linestyle='dotted'
plt.plot(lsxx, muxx, linestyle='', c = 'black',marker='o',ms=4) #linestyle='dotted'
plt.plot(lsxx1, muxx1, linestyle='', c = 'grey',marker='o',ms=4) #linestyle='dotted'
l = np.zeros(20)
mu = np.zeros(20)
lk=20
for i in range(20):
    l[i]=i*10+5
    mu[i]=-0.34+math.log(1/(l[i]+lk))
    #print(p[i],p[i]*l[j])
plt.plot(l, mu,c='grey',label='$l_{k}=20$') #linestyle='dotted'
plt.legend(loc="upper right",fontsize=11,frameon=False)

#plt.vlines(-2.65, 0, 1, linestyles ="dotted", colors ="k")

    
fig.tight_layout()
plt.show()
  
 




