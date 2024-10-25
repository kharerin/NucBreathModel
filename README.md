# NucBreathModel
The codes to generate the data and figures in the paper are provided below:

## Figure 1
The data are stored in the folder yeastsc_lnZ/H4_Input_MNase_200U_Rep1 and were generated using MNase data provided in GEO: ID GSE83123, dataset_H4_Input_MNase_200U Replicate_2. The pipeline for the analyses (extraction and refernce mapping) of the seqeuncing reads are provided in the paper associated with the same GEO ID. To generate the figure 1, run on the terminal:
```
$ python3 Figure_1_mnase200U.py.
```

## Figure 2
To compute density and occupancy run the code hung_soft_lnZ.c found in the folder density_fixedsize_vs_mu. On terminal type the commands:
```
$ gcc hung_soft_lnZ.c -lm 
$ ./a.out.
```
To calculate the critical chemical potential vs wrapping length (the third panel of figure 2) for lk=0 and lk=20, run the codes /critical_mu/hung_soft_lnZ.c and /critical_mu_lk_20/hung_soft_lnZ.c, respectively. To generate the figure 2, run fixedtonks_fig2.py.
