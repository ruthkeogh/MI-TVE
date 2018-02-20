# MI-TVE
Multiple imputation in Cox regression when there are time-varying effects of exposures

The file generate_data.R simulates time-to-event data with time-varying effects of covariates, as described in the paper:
Keogh RH, Morris TP. Multiple imputation in Cox regression when there are time-varying effects of exposures. https://arxiv.org/abs/1706.09187

The file smcfcs_tve provides the code for the MI-TVE-SMC methods described in the paper, and is modified from Jonathan Bartlett's 'smcfcs' R package (https://github.com/jwb133/smcfcs, and also available on CRAN).

The file methods.R provides examples of how to implement the methods described in the above paper by applying them to the simulated data created in generate_data.R. 

