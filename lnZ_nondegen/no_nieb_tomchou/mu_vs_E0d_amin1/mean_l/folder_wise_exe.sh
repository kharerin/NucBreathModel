#!/bin/bash

for i in {201..250}
do
  mkdir nucdir$i
  cp hung_soft_lnZ_set.c ./nucdir$i/hung_soft_lnZ_set.c
  cp jop_script.sh ./nucdir$i/jop_script.sh
  cd nucdir$i
  sed -i "s:replacepath:$(pwd):" hung_soft_lnZ_set.c
  mv hung_soft_lnZ_set.c hung_soft_lnZ_set$i.c
  sed -i "s/tn=1/tn=$i/" hung_soft_lnZ_set$i.c
  sed -i "s/hung_soft_lnZ_set.c/hung_soft_lnZ_set$i.c/" jop_script.sh
  sed -i "s/ms0p1x/coop0_$i/" jop_script.sh
  sed -i "s/den.out/den$i.out/" jop_script.sh
  sed -i "s/denout.txt/denout$i.txt/" jop_script.sh
  ./jop_script.sh
 cd ..
done 
