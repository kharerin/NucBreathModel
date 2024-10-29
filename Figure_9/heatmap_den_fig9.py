import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import string

fig, ax = plt.subplots(3,2,figsize=(7.2,9))

mu2=[]
E02=[]
dg=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_vs_E0d/mean_dg/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    dg.append(float(row[2]))
f.close() 
   
k=0   
for i in range(250):
    for j in range(250):
        if dg[j+k]<91:
           data3[249-i][j]=1169
        else:
           data3[249-i][j]=dg[j+k]    
    k=k+250
    
mu1x=[]
E01x=[]
mu1y=[]
E01y=[]
mu1z=[]
E01z=[]
for i in range(250*250):
    if dg[i]>154 and dg[i]<156:
       mu1x.append(mu2[i])
       E01x.append(E02[i]) 
for i in range(250*250):
    if dg[i]>174 and dg[i]<176:
       mu1y.append(mu2[i])
       E01y.append(E02[i])
for i in range(250*250):
    if dg[i]>164 and dg[i]<166:
       mu1z.append(mu2[i])
       E01z.append(E02[i])
       print('>>>>>>1st model...mu,e0,dg',mu2[i],E02[i],dg[i])

n=-56
E11x=np.multiply(E01x,n)  
E11y=np.multiply(E01y,n)   
E11z=np.multiply(E01z,n)
            
ax1=plt.subplot(3, 2, 1)
ax1.text(-0.12, 1.05, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
#plt.title('NRL',fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_0$',fontsize=12)

plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=20, vmax=400, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
print('min...',data3.min())
plt.scatter(mu1x,E11x, s=3, c = 'k')
plt.scatter(mu1y,E11y, s=3, c = 'grey')
plt.scatter(mu1z,E11z, s=3, c = 'm')
plt.scatter(-12.8,-0.057*56, marker = '+', s=100, c = 'w',zorder=1)

mu2=[]
E02=[]
pO=[]
lm=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
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
plt.plot(mu1x,E11x, linewidth=2, c = 'w',linestyle='dotted')


mu2=[]
pO2=[]
dg2=[]
lm2=[]
data1 = np.zeros( (6, 575) )
data4 = np.zeros( (6, 575) )
data5 = np.zeros( (6, 575) )
i=0
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_vs_E0d/mean_dg_vs_mu/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    pO2.append(float(row[2]))
    dg2.append(float(row[3]))
    lm2.append(float(row[4]))
f.close()
mu1=mu2[0:575]
n=0.0001*20
pO1=np.multiply(pO2,n)
k=0     
for i in range(6):
    for j in range(575):
        data1[5-i][j]=pO1[j+k]
        data4[5-i][j]=dg2[j+k]
        data5[5-i][j]=lm2[j+k]
    k=k+575  

ax1=plt.subplot(3, 2, 2)
ax1.text(-0.15, 1.05, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("densityxa")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('NRL',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,400)
plt.scatter(mu1, data4[0][:], s=3, c = 'brown')
plt.scatter(mu1, data4[1][:], s=3,c = 'r')
plt.scatter(mu1, data4[2][:], s=3,c = 'g')
plt.scatter(mu1, data4[3][:], s=3,c = 'c')
plt.scatter(mu1, data4[4][:], s=3,c = 'b')
plt.scatter(mu1, data4[5][:], s=3,c = 'k')
plt.axhline(y = 91, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 147, color = 'grey', linestyle = '-') 

#panel1

mu2=[]
E02=[]
dg=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d/mean_dg/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    dg.append(float(row[2]))
f.close() 

k=0   
for i in range(250):
    for j in range(250):
        if dg[j+k]<111:
           data3[249-i][j]=1169
        else:
           data3[249-i][j]=dg[j+k]
    k=k+250

mu1x=[]
E01x=[]
mu1y=[]
E01y=[]
mu1z=[]
E01z=[]
for i in range(250*250):
    if dg[i]>154 and dg[i]<156:
       mu1x.append(mu2[i])
       E01x.append(E02[i]) 
for i in range(250*250):
    if dg[i]>174 and dg[i]<176:
       mu1y.append(mu2[i])
       E01y.append(E02[i])
for i in range(250*250):
    if dg[i]>164 and dg[i]<166:
       mu1z.append(mu2[i])
       E01z.append(E02[i])
       print('>>>>>>2nd...mu,e0,dg',mu2[i],E02[i],dg[i]) 

n=-56
E11x=np.multiply(E01x,n)  
E11y=np.multiply(E01y,n)   
E11z=np.multiply(E01z,n)    


ax1=plt.subplot(3, 2, 3)
ax1.text(-0.12, 1.05, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
#plt.title(r'Density ($\rho$)',fontsize=12)
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_0$',fontsize=12)

plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=20, vmax=400, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
print('min...',data3.min())
plt.scatter(mu1x,E11x, s=3, c = 'k')
plt.scatter(mu1y,E11y, s=3, c = 'grey')
plt.scatter(mu1z,E11z, s=3, c = 'm')
plt.scatter(-9.2,-0.064*56, marker = '+', s=100, c = 'w',zorder=1)

mu2=[]
E02=[]
pO=[]
lm=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d/mean_l/d1_mu_E0d.txt','r')
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
plt.plot(mu1x,E11x, linewidth=2, c = 'w',linestyle='dotted')


mu2=[]
pO2=[]
dg2=[]
lm2=[]
data1 = np.zeros( (6, 575) )
data4 = np.zeros( (6, 575) )
data5 = np.zeros( (6, 575) )
i=0
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d/mean_dg_vs_mu/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    pO2.append(float(row[2]))
    dg2.append(float(row[3]))
    lm2.append(float(row[4]))
f.close()
mu1=mu2[0:575]
n=0.0001*20
pO1=np.multiply(pO2,n)
k=0     
for i in range(6):
    for j in range(575):
        data1[5-i][j]=pO1[j+k]
        data4[5-i][j]=dg2[j+k]
        data5[5-i][j]=lm2[j+k]
    k=k+575 

ax1=plt.subplot(3, 2, 4)
ax1.text(-0.15, 1.05, string.ascii_uppercase[3], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("densityxa")
#plt.xlabel('$\mu$',fontsize=12)
plt.ylabel('NRL',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,400)
plt.scatter(mu1, data4[0][:], s=3, c = 'brown')
plt.scatter(mu1, data4[1][:], s=3,c = 'r')
plt.scatter(mu1, data4[2][:], s=3,c = 'g')
plt.scatter(mu1, data4[3][:], s=3,c = 'c')
plt.scatter(mu1, data4[4][:], s=3,c = 'b')
plt.scatter(mu1, data4[5][:], s=3,c = 'k')
plt.axhline(y = 91+20, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 147+20, color = 'grey', linestyle = '-') 

#panel2

mu2=[]
E02=[]
dg=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d_amin1/mean_dg/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    dg.append(float(row[2]))
f.close()
k=0   
for i in range(250):
    for j in range(250):
        if dg[j+k]<20:
           data3[249-i][j]=1169
        else:
           data3[249-i][j]=dg[j+k]
    k=k+250

mu1x=[]
E01x=[]
mu1y=[]
E01y=[]
mu1z=[]
E01z=[]
for i in range(250*250):
    if dg[i]>154 and dg[i]<156:
       mu1x.append(mu2[i])
       E01x.append(E02[i])
for i in range(250*250):
    if dg[i]>174 and dg[i]<176:
       mu1y.append(mu2[i])
       E01y.append(E02[i])
for i in range(250*250):
    if dg[i]>164 and dg[i]<166:
       mu1z.append(mu2[i])
       E01z.append(E02[i])
       print('>>>>>>3rd...mu,e0,dg',mu2[i],E02[i],dg[i])

n=-56
E11x=np.multiply(E01x,n)  
E11y=np.multiply(E01y,n)   
E11z=np.multiply(E01z,n)             


ax1=plt.subplot(3, 2, 5)
ax1.text(-0.12, 1.05, string.ascii_uppercase[4], transform=ax1.transAxes, size=20, weight='bold')
#plt.title(r'Density ($\rho$)',fontsize=12)
plt.xlabel(r'$\mu$',fontsize=12)
plt.ylabel(r'$\epsilon_0$',fontsize=12)

plt.xlim(-60, 40)
#plt.ylim(25,-25)
ax1.set_yticks([-20, -10, 0, 10, 20]) 
ax1.set_yticklabels(['$-0.36$', '$-0.18$', '0','$0.18$', '$0.36$']) 

imx1=plt.imshow(data3,extent=[-60, 40,25,-25], vmin=20, vmax=400, cmap='jet',aspect=2);
plt.colorbar(imx1,shrink=0.5)
print('min...',data3.min())
plt.scatter(mu1x,E11x, s=3, c = 'k',label='NRL=155')
plt.scatter(mu1z,E11z, s=3, c = 'm',label='NRL=165')
plt.scatter(mu1y,E11y, s=3, c = 'grey',label='NRL=175')
plt.scatter(-11.6,-0.096*56, marker = '+', s=100, c = 'w',zorder=1)
plt.legend(loc="upper right",bbox_to_anchor=(1.1, -0.2),borderpad=0.05,labelspacing=0.25,ncol=1,facecolor='w',framealpha=1,frameon=True,fontsize=10) #framealpha=1

mu2=[]
E02=[]
pO=[]
lm=[]
data1 = np.zeros( (250, 250) )
data2 = np.zeros( (250, 250) )
data3 = np.zeros( (250, 250) )
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d_amin1/mean_l/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    E02.append(float(row[1]))
    lm.append(float(row[7]))
f.close() 
j=0
mu1x=[]
E01x=[]
for i in range(250*250):
    if lm[i]>73 and lm[i]<75:
       mu1x.append(mu2[i])
       E01x.append(E02[i])
       j=j+1  
n=-56
E11x=np.multiply(E01x,n)  
plt.plot(mu1x,E11x, linewidth=2, c = 'w',linestyle='dotted')


mu2=[]
pO2=[]
dg2=[]
lm2=[]
data1 = np.zeros( (6, 575) )
data4 = np.zeros( (6, 575) )
data5 = np.zeros( (6, 575) )
i=0
f = open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d_amin1/mean_dg_vs_mu/d1_mu_E0d.txt','r')
for row in f:
    row = row.split(' ')
    mu2.append(float(row[0]))
    pO2.append(float(row[2]))
    dg2.append(float(row[3]))
    lm2.append(float(row[4]))
f.close()
mu1=mu2[0:575]
n=0.0001*20
pO1=np.multiply(pO2,n)
k=0     
for i in range(6):
    for j in range(575):
        data1[5-i][j]=pO1[j+k]
        data4[5-i][j]=dg2[j+k]
        data5[5-i][j]=lm2[j+k]
    k=k+575 

ax1=plt.subplot(3, 2, 6)
ax1.text(-0.15, 1.05, string.ascii_uppercase[5], transform=ax1.transAxes, size=20, weight='bold')
#plt.title("densityxa")
plt.xlabel(r'$\mu$',fontsize=12)
plt.ylabel('NRL',fontsize=12)
plt.xlim(-75, 40)
plt.ylim(0,400)
plt.scatter(mu1, data4[0][:], s=3, c = 'brown',label=r"$\epsilon_{0}=-0.411$")
plt.scatter(mu1, data4[1][:], s=3,c = 'r',label=r"$\epsilon_{0}=-0.274$")
plt.scatter(mu1, data4[2][:], s=3,c = 'g',label=r"$\epsilon_{0}=-0.137$")
plt.scatter(mu1, data4[3][:], s=3,c = 'c',label=r"$\epsilon_{0}=-0.034$")
plt.scatter(mu1, data4[4][:], s=3,c = 'b',label=r"$\epsilon_{0}=0$")
plt.scatter(mu1, data4[5][:], s=3,c = 'k',label=r"$\epsilon_{0}=0.068$")
plt.axhline(y = 1+20, color = 'grey', linestyle = 'dotted') 
plt.axhline(y = 147+20, color = 'grey', linestyle = '-') 
plt.legend(loc="upper right",bbox_to_anchor=(1.05, -0.2),borderpad=0.05,labelspacing=0.25,ncol=2,facecolor='w',framealpha=1,frameon=True,fontsize=10) #framealpha=1
 
fig.tight_layout()
plt.subplots_adjust(wspace=0.22, hspace=0.25)
plt.show()
  
 




