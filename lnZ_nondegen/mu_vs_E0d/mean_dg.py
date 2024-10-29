import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import pandas as pd
import os 

# N is the number of folders corresponding to either E0 or mu is varied N times
N=250
with open('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_vs_E0d/mean_dg/d1_mu_E0d.txt', 'w') as newfile:
     for i in range(N):
         x=i+1
         #fx=''.join((('/nucdir%d'% x),'/d1_mu_E0d.txt'))
         print(x)
         with open(os.path.join(('/media/hungyo/edrive/hung_ens_laptop/hungyo_nucdyn/NucleoModel/lnZ_nondegen/mu_vs_E0d/mean_dg/nucdir%d/'% x), 'd1_mu_E0d.txt')) as f: 
              contents = f.read()
              newfile.write(contents)
  
 




