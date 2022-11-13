% This is the top-level Matlab file to run for the repo model-inference-covid19-nz2021

clear all
close all

% ABCSMC algorithm will be run fitting to data up to the selected 'fittedToDate'
% Vaccine data will be read in from a filename with the datestamp 'vaxDataDate'
% Three options below, corresponding to the three sets of results in the
% aritcle
fittedToDate = datenum('09JAN2022');  vaxDataDate = datenum('17JAN2022');
%fittedToDate = datenum('28NOV2021');  vaxDataDate = datenum('25NOV2021');
%fittedToDate = datenum('24OCT2021');  vaxDataDate = datenum('14OCT2021');

% Hyperparameters for ABCSMC algorithm,
hyperPar.nParticles = 500;      % Number of particles required
hyperPar.acceptRate = 0.5;      % Acceptance rate for each SMC step
hyperPar.eps0 = 6500;           % Starting tolerance (helps early reject candidates with high distance)  6500, 3000, 600 for 9 Jan, 28 Nov, 24 Oct
hyperPar.epsTarg = 1200;        % Target tolerance                                                       1200, 550, 100 for 9 Jan, 28 Nov, 24 Oct
hyperPar.fAcc = 0.005;          % Minimum MCMC acceptance rate
hyperPar.c = 0.01;              % Tuning parameter controlling the minimum proportion of unique particles
hyperPar.S0 = 10;               % Initial number of MCMC trials

% Load data
fName = 'data/epiData.xlsx';          % 2022-01-11.xlsx
fprintf('Loading epi data:  %s\n', fName); 
epiData = readtable(fName);

% Get model parameters
par = getPar(vaxDataDate);                    

% Run ABCSMC inference algorithm if results don't already exist
fName = sprintf('results/ABCSMC_cases_to_%s.mat', datestr(fittedToDate));
if exist(fName) == 0
   fitToData(epiData, fittedToDate, par, hyperPar);
end

% Run forward model if results don't already exist
fName = sprintf('results/results_fit_to_%s.mat', datestr(fittedToDate));
if exist(fName) == 0
   runForwardSims(epiData, fittedToDate, par);
end

% Plot graphs
plotGraphs(epiData, fittedToDate);


