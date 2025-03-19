README for BTSPvsSTDP, created 12/2024

This code is associated with the article " Synaptic plasticity rules driving representational shifting in the hippocampus ", Nature Neuroscience (2025) by Antoine Madar, Anqi Jiang, Can Dong and Mark Sheffield (University of Chicago).

Contact: madar@uchicago.edu or antoinemadar461@gmail.com

Figures 1a-c, 5, 6, 7a-c and associated supplementary figures: analysis of experimental data from Dong et al. 2021 (https://www.nature.com/articles/s41467-021-23260-3)

The data loaded in the following Matlab scripts is available at COMtraj/Main/data/CA1_COM_withPFs.mat and CA3_COM_withPFs.mat.
The main routine is: COMtraj/Main/CanDataAnalysis.m
In addition to plots, it can create some matfiles, such as CA1_SlopeDistrib.mat and CA1perAnimal_slopes_props.mat, that are reused in other scripts.
For visualizing example PFs and additional plots: COMtraj/Main/CanDataAnalysis_ExPFs.m
Stats for Fig 1b are done in: COMtraj/Main/CanDataAnalysis_LinRegBootstrapStats.m
Stats for Fig 6f-g and S19b are done in: COMtraj/Main/CanDataAnalysis_NonLinRegBootstrapStats.m

Figures 1d: Optogenetic experiment

Data tables loaded in the following routines are in COMtraj/Opto/data
To analyze PF emerged in sessions when opto silencing was performed at beginning of session (opto_first): COMtraj/Opto/OptoCOMtraj.m 
To analyze PF emerged in sessions when opto silencing was performed in the middle of session: COMtraj/Opto/OptoCOMtraj_OptoLater_PFinOpto.m 
To pool and analyze both sets of PFs: COMtraj/Opto/OptoCOMtraj_Pooled.m

Figures 2 and associated supplementary figures: STDP models

-	Model core: STDP/STDPplus_LIFadapt.m 
-	Fig 2c-d, S5-6: STDP/STDPmodelLIF_Repeat1Rule_BATCH.m
-	Fig 2e: STDP/SummaryModels_vs_CA1data.m
-	Fig 2f, S3 and S7-10: STDP/STDPmodelLIF_BATCH_ParamSpaceRepeat.m
-	Fig S11a-b: STDP/triplet-FigS11/InVitro_STDPtripletProtocol.m 
-	Fig S11c-h: STDP/ triplet-FigS11/STDPtripletLIF_Repeat1Rule_BATCH.m (which uses the function STDPtriplet.m)

Figures 3, 4, 7c and associated supplementary figures: heterosynaptic BTSP plasticity models

-	Model core: BTSP/BTSPplus_LIFadapt.m (Params.pCSdynamic = 0)
-	Fig 3b-f: BTSP/BTSPmodelLIF_Repeat1Rule_BATCH.m (Params.pCSdynamic = 0)
-	Fig 3g-h: BTSP/SummaryBTSPmodels_vs_CA1data.m
-	Fig 4: BTSP/BTSPmodelLIF_Repeat1Rule_BATCH.m + SummaryBTSPmodels_vs_CA3data.m
-	Fig S12: BTSP/InVitroBittner
-	Fig S13: BTSP/MilsteinExp/InVivoMilstein2021_1CS.m
-	Fig S14-15: BTSP/MilsteinExp/InVivoMilstein2021_1CSbatch.m
-	Fig S14h-j and 15h-i: BTSP/MilsteinExp/InVivoMilstein2021_1CSbatch_PampSpace.m
-	Fig 7c bootstrap stats and Fig S20: SummaryBTSPmodels_vs_CA1data_MSD.m

Fig S16: BTSP as a weight-dependent homosynaptic plasticity rule

-	Model core: BTSP/biBTSP/MilsteinBTSP_LIFadapt_ODEs.m or BTSP_MilsteinModel_LIFadapt.m (same as the ODE-only function but includes input simulation)
-	Fig S16a-c: BTSP/biBTSP/InVivoMilstein2021_MilsteinModel.m
-	Fig S16d-e: BTSP/biBTSP/InVivoMilstein2021_MilsteinModel_batch.m
Both scripts above use the function MilsteinBTSP_LIFadapt_ODEs.m
-	Simulations for Fig S16f-i: BTSP/biBTSP/BTSP_MilsteinModel_LIF_Repeat1Rule_BATCH.m
This last script uses the function BTSP_MilsteinModel_LIFadapt.m
-	Fig S16g-i: BTSP/SummaryBTSPmodels_vs_CA1data.m

Figure 7d and S21: BTSP models with dynamic p(CS)

-	Model core: BTSP/BTSPplus_LIFadapt.m (Params.pCSdynamic = 1)
-	Simulations: BTSP/BTSPmodelLIF_Repeat1Rule_BATCH.m (Params.pCSdynamic = 1)
-	Figures: BTSP/SummaryBTSPmodels_vs_CA1data.m

*Dependencies*

- Brewermap : https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer-attractive-and-distinctive-colormaps
- plot_ci : https://www.mathworks.com/matlabcentral/fileexchange/31752-plot_ci
- violinplot: Bechtold, Bastian, 2016. https://github.com/bastibe/Violinplot-Matlab, DOI: 10.5281/zenodo.4559847
