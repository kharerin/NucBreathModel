# NucBreathModel
The codes to generate the data and figures in the paper are provided below:

## Figure 1
Go to folder Figure_1
```
$ cd Figure_1
```
To generate the figure 1, run on the terminal:
```
$ python3 Figure_1_mnase200U.py
```
The data is stored in the folder yeastsc_lnZ/H4_Input_MNase_200U_Rep1 and were generated using MNase data provided in GEO: ID GSE83123, dataset_H4_Input_MNase_200U Replicate_2. The pipeline for the analyses (extraction and refernce mapping) of the seqeuncing reads are provided in the paper associated with the same GEO ID.

## Figure 2
To generate the figure 2, run:
```
$ python3 fixedtonks_fig2.py
```
Data has been pregenerated for you and can be found in the same folder. Otherwise it can be generated using the commands:
```
$ gcc hung_soft_lnZ.c -lm 
$ ./a.out
```
It computs density and occupancy and be found in the folder density_fixedsize_vs_mu. To calculate the critical chemical potential vs wrapping length (the third panel of figure 2) for $l_{k}=0$ and $l_k=20$, run the codes /critical_mu/hung_soft_lnZ.c and /critical_mu_lk_20/hung_soft_lnZ.c, respectively. 

## Figure 3
To generate the figure 3, run:
```
$ python3 heatmap_den_fig3.py
```
Data has been pregenerated for you and can be found in the same folder. Otherwise it can be generated using the following commands.
```
$ cd lnZ_nondegen/mu_vs_E0d/mean_l/
$ gcc hung_soft_lnz.c -lm
$ ./a.out
```
Come back to initial directory and run:
```
$ cd lnZ_nondegen/no_nieb_cedric/
$ gcc equi_density_cedric.c -lm
$ ./a.out
```
Come back to initial directory and run:
```
$ cd lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d/mean_l/
$ gcc hung_soft_lnZ.c -lm
$ ./a.out
```
Come back to initial directory and run:
```
$ cd lnZ_nondegen/no_nieb_tomchou/mu_vs_E0d_amin1/mean_l/
$ ./folder_wise_exe.sh
$ cd ..
$ python3 mean_l.py
```
Come back to initial directory and run:
```
$ cd lnZ_nondegen/fixed_tonks_gas_model/mu_vs_E0d/
$ gcc hung_soft_lnz.c -lm
$ ./a.out
```

