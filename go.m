clear all
close all

% ABCSMC algorithm will be run using data up to this date:
% Vaccine data will be read in from a filename with this datestamp
fittedToDate = datenum('09JAN2022');  vaxDataDate = datenum('17JAN2022');
%fittedToDate = datenum('28NOV2021');  vaxDataDate = datenum('25NOV2021');
%fittedToDate = datenum('24OCT2021');  vaxDataDate = datenum('14OCT2021');


% Load data
fName = 'data/epiData.xlsx';          % 2022-01-11.xlsx
fprintf('Loading epi data:  %s\n', fName); 
epiData = readtable(fName);

% Get model parameters
par = getPar(vaxDataDate);                    

% Run ABCSMC inference algorithm if results don't already exist
fName = sprintf('results/ABCSMC_cases_to_%s.mat', datestr(fittedToDate));
if exist(fName) == 0
   runABC(epiData, fittedToDate, par);
end

% Run forward model if results don't already exist
fName = sprintf('results/results_fit_to_%s.mat', datestr(fittedToDate));
if exist(fName) == 0
   runForwardModel(epiData, fittedToDate, par);
end

% Plot graphs
plotGraphs(epiData, fittedToDate);


