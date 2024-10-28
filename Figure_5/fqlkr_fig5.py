import matplotlib.pyplot as plt
import numpy as np
import math
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import pandas as pd
import string
    
fig, ax = plt.subplots(3,2,figsize=(7,10))

f1='/lnZ_nondegen/no_nieb_cedric/' # model 1

ax2=plt.subplot(3, 2, 1)
ax2.text(-0.15, 1.05, string.ascii_uppercase[0], transform=ax2.transAxes, size=20, weight='bold')
axins = ax2.inset_axes([0.3, 0.55, 0.36, 0.37])
#axins.set_title('log(y) vs x')
axins.set_yscale('log')
plt.title(r"$\epsilon_0=-0.36$")
#plt.xlabel('Fragments (bp)',fontsize=11)
plt.ylabel('Frequency',fontsize=11)
plt.axis([91, 147, 0, 0.05])
#plt.yscale("log")

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   
    
f2='E0d_-0p36kT_mu_-25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'black',label='')
axins.plot(o11, o3, c = 'black')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    

f2='E0d_-0p36kT_mu_-15kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkolivegreen',label='')
axins.plot(o11, o3, c = 'darkolivegreen')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i 

f2='E0d_-0p36kT_mu_-6kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'salmon',linestyle='dashed',label='')
axins.plot(o11, o3, c = 'salmon',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   

f2='E0d_-0p36kT_mu_-5p25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'b',linestyle='dashed',label=r"$\mu*=-5.25$")
axins.plot(o11, o3, c = 'b',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   

f2='E0d_-0p36kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkred',label='')
axins.plot(o11, o3, c = 'darkred')
#plt.setp(axins.get_xticklabels(), visible=False)
#plt.setp(axins.get_yticklabels(), visible=False)
#axins.set_xticks([])
#axins.set_yticks([])
plt.legend(loc="lower right",borderpad=0.25,labelspacing=0.25,fontsize=10) #bbox_to_anchor=(1.03, -0.04)
#axins.set_xlim(91, 147)
#axins.set_ylim(0, 0.05)
#plt.setp(ax2.get_xticklabels(), visible=False)


ax2=plt.subplot(3, 2, 2)
ax2.text(-0.15, 1.05, string.ascii_uppercase[1], transform=ax2.transAxes, size=20, weight='bold')
axins = ax2.inset_axes([0.3, 0.55, 0.36, 0.37])
#axins.set_title('log(y) vs x')
axins.set_yscale('log')
plt.title(r"$\epsilon_0=-0.18$")
#plt.xlabel('Fragments (bp)',fontsize=11)
#plt.ylabel('Frequency',fontsize=11)
plt.axis([91, 147, 0, 0.05])
#plt.yscale("log")

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i 

f2='E0d_-0p18kT_mu_-25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')   
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'black',label='')
axins.plot(o11, o3, c = 'black')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i  

f2='E0d_-0p18kT_mu_-15kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx) 
f = open(fx,'r')   
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkolivegreen',label='')
axins.plot(o11, o3, c = 'darkolivegreen')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   

f2='E0d_-0p18kT_mu_-6kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'salmon',linestyle='dashed',label='')
axins.plot(o11, o3, c = 'salmon',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i

f2='E0d_-0p18kT_mu_-5p76kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx) 
f = open(fx,'r')     
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'b',linestyle='dashed',label=r"$\mu*=-5.76$")
axins.plot(o11, o3, c = 'b',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   

f2='E0d_-0p18kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)   
f = open(fx,'r')   
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkred',label='')
axins.plot(o11, o3, c = 'darkred')
#plt.setp(axins.get_xticklabels(), visible=False)
#plt.setp(axins.get_yticklabels(), visible=False)
#axins.set_xticks([])
#axins.set_yticks([])
plt.legend(loc="lower right",borderpad=0.25,labelspacing=0.25,fontsize=10) #bbox_to_anchor=(1.03, -0.04)

#axins.set_xlim(91, 147)
#axins.set_ylim(0, 0.05)
#plt.setp(ax2.get_xticklabels(), visible=False)


ax2=plt.subplot(3, 2, 3)
ax2.text(-0.15, 1.05, string.ascii_uppercase[2], transform=ax2.transAxes, size=20, weight='bold')
axins = ax2.inset_axes([0.3, 0.55, 0.36, 0.37])
#axins.set_title('log(y) vs x')
axins.set_yscale('log')
plt.title(r"$\epsilon_0=-0.09$")
#plt.xlabel('Fragments (bp)',fontsize=11)
plt.ylabel('Frequency',fontsize=11)
plt.axis([91, 147, 0, 0.05])
#plt.yscale("log")

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i 

f2='E0d_-0p09kT_mu_-25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')     
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'black',label='')
axins.plot(o11, o3, c = 'black')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i  

f2='E0d_-0p09kT_mu_-15kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')      
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkolivegreen',label='')
axins.plot(o11, o3, c = 'darkolivegreen')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i 

f2='E0d_-0p09kT_mu_-6kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')      
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'salmon',linestyle='dashed',label='')
axins.plot(o11, o3, c = 'salmon',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    
    
f2='E0d_-0p09kT_mu_-6p46kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')  
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'b',linestyle='dashed',label=r"$\mu*=-6.46$")
axins.plot(o11, o3, c = 'b',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    
    
f2='E0d_-0p09kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')  
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkred',label='')
axins.plot(o11, o3, c = 'darkred')
#plt.setp(axins.get_xticklabels(), visible=False)
#plt.setp(axins.get_yticklabels(), visible=False)
#axins.set_xticks([])
#axins.set_yticks([])
plt.legend(loc="lower right",borderpad=0.25,labelspacing=0.25,fontsize=10) #bbox_to_anchor=(1.03, -0.04)
#axins.set_xlim(91, 147)
#axins.set_ylim(0, 0.05)
#plt.setp(ax2.get_xticklabels(), visible=False)


ax2=plt.subplot(3, 2, 4)
ax2.text(-0.15, 1.05, string.ascii_uppercase[3], transform=ax2.transAxes, size=20, weight='bold')
axins = ax2.inset_axes([0.3, 0.55, 0.36, 0.37])
#axins.set_title('log(y) vs x')
axins.set_yscale('log')
plt.title(r"$\epsilon_0=0$")
#plt.xlabel('Fragments (bp)',fontsize=11)
#plt.ylabel('Frequency',fontsize=11)
plt.axis([91, 147, 0, 0.05])
#plt.yscale("log")

#o1=[]; o2=[]; o3=[]; sum1=0;
#o11 = np.zeros(57)
#for i in range(57):
#    o11[i]=91+i    
#f2='E0d_-0kT_mu_-25kT/frag_nie2_occ_0pxx.txt'
#fx=''.join((f1,f2))
#print('filename...',fx)  
#f = open(fx,'r')
#for row in f:
#    row = row.split(' ')
#    o1.append(float(row[0]))
#    o2.append(float(row[1]))
#    o3.append(float(row[2]))
#f.close()
#plt.plot(o11, o3, c = 'purple',label="$\mu$=-25")
#axins.plot(o11, o3, c = 'purple')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i  

f2='E0d_-0kT_mu_-15kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')      
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkolivegreen',label='')
axins.plot(o11, o3, c = 'darkolivegreen')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i  

f2='E0d_-0kT_mu_-6kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'salmon',linestyle='dashed',label='')
axins.plot(o11, o3, c = 'salmon',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    
    
f2='E0d_-0kT_mu_-10p4kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r') 
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'b',linestyle='dashed',label=r"$\mu*=-10.4$")
axins.plot(o11, o3, c = 'b',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    
    
f2='E0d_-0kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r') 
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkred',label='')
axins.plot(o11, o3, c = 'darkred')
#plt.setp(axins.get_xticklabels(), visible=False)
#plt.setp(axins.get_yticklabels(), visible=False)
#axins.set_xticks([])
#axins.set_yticks([])
plt.legend(loc="lower right",borderpad=0.25,labelspacing=0.25,fontsize=10) #bbox_to_anchor=(1.03, -0.04)
#axins.set_xlim(91, 147)
#axins.set_ylim(0, 0.05)
#plt.setp(ax2.get_xticklabels(), visible=False)

ax2=plt.subplot(3, 2, 5)
ax2.text(-0.15, 1.05, string.ascii_uppercase[4], transform=ax2.transAxes, size=20, weight='bold')
axins = ax2.inset_axes([0.3, 0.55, 0.36, 0.37])
#axins.set_title('log(y) vs x')
axins.set_yscale('log')
plt.title(r"$\epsilon_0=0.09$")
plt.xlabel('Fragments (bp)',fontsize=11)
plt.ylabel('Frequency',fontsize=11)
plt.axis([91, 147, 0, 0.05])
#plt.yscale("log")

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    

f2='E0d_0p09kT_mu_-25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'black',label=r"$\mu=-25$")
axins.plot(o11, o3, c = 'black')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i  

f2='E0d_0p09kT_mu_-15kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkolivegreen',label=r"$\mu=-15$")
axins.plot(o11, o3, c = 'darkolivegreen')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i 

f2='E0d_0p09kT_mu_-6kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'salmon',linestyle='dashed',label=r"$\mu=-6$")
axins.plot(o11, o3, c = 'salmon',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   

f2='E0d_0p09kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'darkred',label=r"$\mu=25$")
axins.plot(o11, o3, c = 'darkred')
#plt.setp(axins.get_xticklabels(), visible=False)
#plt.setp(axins.get_yticklabels(), visible=False)
#axins.set_xticks([])
#axins.set_yticks([])
plt.legend(loc="lower right",bbox_to_anchor=(0.9, -0.5),borderpad=0.25,labelspacing=0.25,ncol=2,fontsize=10) #bbox_to_anchor=(1.03, -0.04)
#axins.set_xlim(91, 147)
#axins.set_ylim(0, 0.05)
#plt.setp(ax2.get_xticklabels(), visible=False)


ax2=plt.subplot(3, 2, 6)
ax2.text(-0.15, 1.05, string.ascii_uppercase[5], transform=ax2.transAxes, size=20, weight='bold')
axins = ax2.inset_axes([0.3, 0.55, 0.36, 0.37])
#axins.set_title('log(y) vs x')
axins.set_yscale('log')
plt.title(r"$\mu=25$")
plt.xlabel('Fragments (bp)',fontsize=11)
#plt.ylabel('Frequency',fontsize=11)
plt.axis([91, 147, 0, 0.05])
#plt.yscale("log")

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    
    
f2='E0d_-0p36kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r') 
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'lime',label=r"$\epsilon_0=-0.36$")
axins.plot(o11, o3, c = 'lime')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i    
    
f2='E0d_-0p18kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r') 
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'cyan',label=r"$\epsilon_0=-0.18$")
axins.plot(o11, o3, c = 'cyan')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i   

f2='E0d_-0p09kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')     
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'lightgreen',linestyle='dashed',label=r"$\epsilon_0=-0.09$")
axins.plot(o11, o3, c = 'lightgreen',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i  

f2='E0d_-0kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')        
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'olive',linestyle='dashed',label=r"$\epsilon_0=0$")
axins.plot(o11, o3, c = 'olive',linestyle='dashed')

o1=[]; o2=[]; o3=[]; sum1=0;
o11 = np.zeros(57)
for i in range(57):
    o11[i]=91+i 

f2='E0d_0p09kT_mu_25kT/frag_nie2_occ_0pxx.txt'
fx=''.join((f1,f2))
print('filename...',fx)  
f = open(fx,'r')    
for row in f:
    row = row.split(' ')
    o1.append(float(row[0]))
    o2.append(float(row[1]))
    o3.append(float(row[2]))
f.close()
plt.plot(o11, o3, c = 'teal',label=r"$\epsilon_0=0.09$")
axins.plot(o11, o3, c = 'teal')
#plt.setp(axins.get_xticklabels(), visible=False)
#plt.setp(axins.get_yticklabels(), visible=False)
#axins.set_xticks([])
#axins.set_yticks([])
plt.legend(loc="lower right",bbox_to_anchor=(0.9, -0.55),borderpad=0.25,labelspacing=0.25,ncol=2,fontsize=10) #bbox_to_anchor=(1.03, -0.04)

#axins.set_xlim(91, 147)
#axins.set_ylim(0, 0.05)
#plt.setp(ax2.get_xticklabels(), visible=False)


fig.tight_layout()
plt.subplots_adjust(wspace=0.2, hspace=0.25)
plt.show()
 




