# model-inference-covid19-nz2021
**Using mechanistic model-based inference to understand and project epidemic dynamics with time-varying contact and vaccination rates: New Zealand’s 2021 Covid-19 outbreak**

Matlab code to reproduce the analysis in "Using mechanistic model-based inference to understand and project epidemic dynamics with time-varying contact and vaccination rates: New Zealand’s 2021 Covid-19 outbreak"

- Run go.m with a specified choice of fittedToDate to fit the model using case data up to that date and plot the main figures in the article. go.m is the top-level Matlab script that calls three functions:
	1. fitToData.m         [if the relevant results file /results/ABCSMC_cases_to_DD-MMM-YYYY.mat does not exist] - runs the ABCSMC algorithm described in the paper to generate a sample from the posterior distribution of the values of the time-varying control function C(t) by fitting to data on new daily cases up to the specified time point
	2. runForwardSims.m [if the relevant results file /results/results_fit_to_DD-MMM-YYYY.mat does not exist]  - takes the posterior samples saved in ABCSMC_cases_to_DD-MMM-YYYY does not exist and runs the forward model for a longer period of time to generate projections of the future epidemic dynamics under the current estmiated values of the control function
	3. plotGraphs.m                                                                                        - takes the results of the forward simulations saved in /results/results_fit_to_DD-MMM-YYYY and plots graphs with data superimposed
- Run drawVaccineGraphs.m to reproduce Fig. 1 in the article.
- Run plotPriorDraws.m to reproduce Supplementary Figure S2.



**RESULTS**

Results are saved in /results/

- results/ABCSMC_cases_to_DD-MMM-YYYY.mat contains the results of the ABCSMC parameter inference part of the method fitted to - data up TO DD-MMM-YYYY.
- results/results_fit_to_DD-MMM-YYYY.mat contains the output of the forward model using the parameters saved in ABCSMC_cases_to_DD-MMM-YYYY 
- figures/results_fit_to_DD-MMM-YYYY.png contains grpahs of the model output with data superimposed



**DATA**

Data files are read in from /data/



**POPULATION DATA**

- data/nzpopdist.xlsx        New Zealand population size in 5-year age bands
- data/nzcontmatrix.xlsx     Contact matrix for the New Zealand popultion from Prem et al.



**EPIDEMIC DATA**

- data/epiData.xlsx          Number of new reported cases of Covid-19 and number of hopsital beds occupied with Covid-19 cases on each day from 17 Aug 2021 to 23 Jan 2022




**VACCINE DATA**

- data/dates_vaccine.csv	                Dates covered by the vaccine data
- data/firstDoseProp_akl_YYYY-MM-DD.csv      Proportion of the population who had received their 1st (but not 2nd) dose in 5-year age bands and on each date shown in dates_vaccine.csv, according to MOH data available on YYYY-MM-DD
- data/secondDoseProp_akl_YYYY-MM-DD.csv     Proportion of the population who had received their 2nd dose in 5-year age bands and on each date shown in dates_vaccine.csv, according to MOH data available on YYYY-MM-DD



**REFF ESTIMATES FROM EPINOW2**

- data/EpiNow2_estimates_AllRegions_MOHdata_YYYYMMDD.csv         Output from the Epinow2 model including estimates of the time-varying reproduciton number 



