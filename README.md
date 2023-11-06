# BTSPvsSTDP
This code is associated with the article "BTSP, not STDP, Drives Shifts in Hippocampal Representations During Familiarization", by Antoine Madar, Can Dong and Mark Sheffield (University of Chicago) available on BioRxiv: https://www.biorxiv.org/content/10.1101/2023.10.17.562791v1

Contact: madar@uchicago.edu or antoinemadar461@gmail.com

*Figures 1, 5-7 and associated supplementary figures *
(analysis of experimental data from https://www.nature.com/articles/s41467-021-23260-3):

-	COMtraj/CanDataAnalysis.m
-	COMtraj/CanDataAnalysis_ExPFs.m

The data analyzed is available in that folder as well (CA1_COM_withPFs.mat and CA3_COM_withPFs.mat)

*Figures 2 and associated supplementary figures* 
(STDP models):

-	Model core: STDP/STDPplus_LIFadapt.m 
-	Fig 2C-D, S3-4: STDP/STDPmodelLIF_Repeat1Rule_BATCH.m
-	Fig 2E: STDP/SummaryModels_vs_CA1data.m
-	Fig 2F, S5-7: STDP/STDPmodelLIF_BATCH_ParamSpaceRepeat.m

*Figures 3, 4 and associated supplementary figures* 
(BTSP models):

-	Model core: BTSP/BTSPplus_LIFadapt.m
-	Fig 3D-H: BTSP/BTSPmodelLIF_Repeat1Rule_BATCH.m
-	Fig 3I-J: SummaryBTSPmodels_vs_CA1data.m
-	Fig 4: BTSP/BTSPmodelLIF_Repeat1Rule_BATCH.m + SummaryBTSPmodels_vs_CA3data.m
-	Fig S8: BTSP/ InVitroBittner
-	Fig S9: BTSP/MilsteinExp/InVivoMilstein2021_1CS.m
-	Fig S10-11: BTSP/MilsteinExp/InVivoMilstein2021_1CSbatch.m
-	Fig S10H-J and Fig 3B-C: BTSP/MilsteinExp/InVivoMilstein2021_1CSbatch_PampSpace.m

*Dependencies*
Brewermap : https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps  
plot_ci : https://www.mathworks.com/matlabcentral/fileexchange/31752-plot_ci
violinplot: Bechtold, Bastian, 2016. https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847
