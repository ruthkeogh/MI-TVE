# MI-TVE
Multiple imputation in Cox regression when there are time-varying effects of exposures

The file generate_data.R simulates time-to-event data with time-varying effects of covariates, as described in the paper:
Keogh RH, Morris TP. Multiple imputation in Cox regression when there are time-varying effects of exposures. https://arxiv.org/abs/1706.09187

The file smcfcs_tve provides the code for the MI-TVE-SMC methods described in the paper, and is modified from Jonathan Bartlett's 'smcfcs' R package (https://github.com/jwb133/smcfcs, and also available on CRAN).

The file methods.R provides examples of how to implement the methods described in the above paper by applying them to the simulated data created in generate_data.R. 

Example analyses of data from the Rotterdam breast cancer data (please refer to our manuscript, section 6) are given in the files rotterdam_example_completedata.Rmd and rotterdam_example_MI_TVE_SMC.Rmd. The Rotterdam data is freely available from https://www.imbi.uni-freiburg.de/Royston-Sauerbrei-book/index.html#datasets. 

The files waldtype_test_mi.R and waldtype_test_mi_rotterdam.R perform the test of Li et al (J Am Stata Assoc 1991; 86: 1065-1073) and Meng & Rubin (Biometrika 1992; 79: 103-111) and are referred to in the R analysis files. 

The files model_selection_procedure.R and model_selection_procedure_MISMC.R perform the MI-MTVE algorithm described in the manuscript. This code is called from the Rotterdam example analysis code. 

