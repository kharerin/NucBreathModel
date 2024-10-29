import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import string
    
fig, ax = plt.subplots(3,2,figsize=(7,10))

f1='/model_1/g_r_lnZ/E0d_e0_amax_amin_hc_a_0/' # data and codes
ax1=plt.subplot(3, 2, 1)
ax1.text(-0.12, 1.05, string.ascii_uppercase[0], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\epsilon_0=-0.36$")
#plt.xlabel('Distance s (kbp)',fontsize=11)
plt.ylabel('$g^{(2)}$(s)',fontsize=11)
plt.xlim(0, 0.5)
plt.ylim(0,10)

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f2='E0d_-20kT_mu_-25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
for i in range(2000):    
    o11[i]=(o1[i])/1000

plt.plot(o11[0:1000], o2[0:1000], c = 'black',label='')

o1=[]; o2=[];
f2='E0d_-20kT_mu_-15kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label='')

o1=[]; o2=[];
f2='E0d_-20kT_mu_-10kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'salmon',linestyle='-',label='')

o1=[]; o2=[];
f2='E0d_-20kT_mu_-5p25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'b',linestyle='dashed',label=r"$\mu*=-5.25$")

o1=[]; o2=[];
f2='E0d_-20kT_mu_0kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'r',label='')

o1=[]; o2=[];
f2='E0d_-20kT_mu_5kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label='')

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="upper right",fontsize=10) #bbox_to_anchor=(1.03, 1.04)


ax1=plt.subplot(3, 2, 2)
ax1.text(-0.12, 1.05, string.ascii_uppercase[1], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\epsilon_0=-0.18$")
#plt.xlabel('Distance s (kbp)',fontsize=11)
plt.ylabel('$g^{(2)}$(s)',fontsize=11)
plt.xlim(0, 0.5)
plt.ylim(0,10)

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f2='E0d_-10kT_mu_-25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
for i in range(2000):    
    o11[i]=(o1[i])/1000

plt.plot(o11[0:1000], o2[0:1000], c = 'black',label='')

o1=[]; o2=[];
f2='E0d_-10kT_mu_-15kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label='')

o1=[]; o2=[];
f2='E0d_-10kT_mu_-10kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'salmon',linestyle='-',label='')

o1=[]; o2=[];
f2='E0d_-10kT_mu_-5p76kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'b',linestyle='dashed',label=r"$\mu*=-5.76$")

o1=[]; o2=[];
f2='E0d_-10kT_mu_0kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'r',label='')

o1=[]; o2=[];
f2='E0d_-10kT_mu_5kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label='')

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="upper right",fontsize=10)


ax1=plt.subplot(3, 2, 3)
ax1.text(-0.12, 1.05, string.ascii_uppercase[2], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\epsilon_0=-0.09$")
#plt.xlabel('Distance s (kbp)',fontsize=11)
plt.ylabel('$g^{(2)}$(s)',fontsize=11)
plt.xlim(0, 0.5)
plt.ylim(0,10)

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f2='E0d_-5kT_mu_-25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
for i in range(2000):    
    o11[i]=(o1[i])/1000

plt.plot(o11[0:1000], o2[0:1000], c = 'black',label='')

o1=[]; o2=[];
f2='E0d_-5kT_mu_-15kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label='')

o1=[]; o2=[];
f2='E0d_-5kT_mu_-10kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'salmon',linestyle='-',label='')

o1=[]; o2=[];
f2='E0d_-5kT_mu_-6p46kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'b',linestyle='dashed',label=r"$\mu*=-6.46$")

o1=[]; o2=[];
f2='E0d_-5kT_mu_0kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'r',label='')

o1=[]; o2=[];
f2='E0d_-5kT_mu_5kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label='')

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="upper right",fontsize=10)


ax1=plt.subplot(3, 2, 4)
ax1.text(-0.12, 1.05, string.ascii_uppercase[3], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\epsilon_0=0$")
#plt.xlabel('Distance s (kbp)',fontsize=11)
plt.ylabel('$g^{(2)}$(s)',fontsize=11)
plt.xlim(0, 0.5)
plt.ylim(0,10)

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f2='E0d_-0kT_mu_-25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
for i in range(2000):    
    o11[i]=(o1[i])/1000

plt.plot(o11[0:1000], o2[0:1000], c = 'black',label='')

o1=[]; o2=[];
f2='E0d_-0kT_mu_-15kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label='')

o1=[]; o2=[];
f2='E0d_-0kT_mu_-10kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'salmon',linestyle='-',label='')

o1=[]; o2=[];
f2='E0d_-0kT_mu_-10p4kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'b',linestyle='dashed',label=r"$\mu*=-10.4$")

o1=[]; o2=[];
f2='E0d_-0kT_mu_0kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()
plt.plot(o11[0:1000], o2[0:1000], c = 'r',label='')

o1=[]; o2=[];
f2='E0d_-0kT_mu_5kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label='')

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="upper right",fontsize=10)
 

ax1=plt.subplot(3, 2, 5)
ax1.text(-0.12, 1.05, string.ascii_uppercase[4], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\epsilon_0=0.09$")
plt.xlabel('Distance s (kbp)',fontsize=12)
plt.ylabel('$g^{(2)}$(s)',fontsize=12)
plt.xlim(0, 0.5)
plt.ylim(0,10)

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f2='E0d_5kT_mu_-25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

for i in range(2000):    
    o11[i]=(o1[i])/1000

plt.plot(o11[0:1000], o2[0:1000], c = 'black',label=r"$\mu=-25$")

o1=[]; o2=[];
f2='E0d_5kT_mu_-15kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'darkolivegreen',label=r"$\mu=-15$")

o1=[]; o2=[];
f2='E0d_5kT_mu_-10kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'salmon',linestyle='-',label=r"$\mu=-10$")

o1=[]; o2=[];
f2='E0d_5kT_mu_0kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'r',label=r"$\mu=0$")

o1=[]; o2=[];
f2='E0d_5kT_mu_5kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'darkred',label=r"$\mu=5$")

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="lower right",bbox_to_anchor=(0.9, -0.55),borderpad=0.25,labelspacing=0.25,ncol=2,fontsize=10) #bbox_to_anchor=(1.03, -0.04)


ax1=plt.subplot(3, 2, 6)
ax1.text(-0.12, 1.05, string.ascii_uppercase[5], transform=ax1.transAxes, size=20, weight='bold')
plt.title(r"$\mu=25$")
plt.xlabel('Distance s (kbp)',fontsize=12)
plt.ylabel('$g^{(2)}$(s)',fontsize=12)
plt.xlim(0, 0.5)
plt.ylim(0,10)

o1=[]; o2=[]; o3=[]; o4=[]; o11=[]
o11 = np.zeros(2000)
f2='E0d_-20kT_mu_25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

for i in range(2000):    
    o11[i]=(o1[i])/1000

plt.plot(o11[0:1000], o2[0:1000], c = 'lime',label=r"$\epsilon_0=-0.36$")

o1=[]; o2=[];
f2='E0d_-10kT_mu_25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'cyan',label=r"$\epsilon_0=-0.18$")

o1=[]; o2=[];
f2='E0d_-5kT_mu_25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'lightgreen',linestyle='dashed',label=r"$\epsilon_0=-0.09$")


o1=[]; o2=[];
f2='E0d_-0kT_mu_25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'olive',linestyle='dashed',label=r"$\epsilon_0=0$")

o1=[]; o2=[];
f2='E0d_5kT_mu_25kT/gr.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
f.close()

plt.plot(o11[0:1000], o2[0:1000], c = 'teal',label=r"$\epsilon_0=0.09$")

#plt.plot(o11[955:1500], o3[955:1500], c = 'blue')
#plt.plot(o11[955:1500], o4[955:1500], c = 'k') 
plt.xticks(np.arange(0, 0.5, 0.1))
plt.yticks(np.arange(0, 10, 2))
plt.legend(loc="lower right",bbox_to_anchor=(0.9, -0.55),borderpad=0.25,labelspacing=0.25,ncol=2,fontsize=10) #bbox_to_anchor=(1.03, -0.04)

fig.tight_layout()
plt.show()
  
 




