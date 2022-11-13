function runForwardModel(epiData, fittedToDate, par)


fNameABC = sprintf('results/ABCSMC_cases_to_%s.mat', datestr(fittedToDate));     
fprintf('   Loading ABC results:    %s\n', fNameABC);
load(fNameABC, 'Cvalues' , 'modelOutput', 'distances');

par.tEnd = 180; % 270 

nReps = 2500;                             % number of repeat simulations to perform
filterFactor = 0.2;                       % retain this proportion of trajectories that fit recent case numbers most closely
filterPeriod = 1+fittedToDate-datenum(epiData.Date(1));     % fit over this period of time (days) up to 'fittedToDate'

ReffQuantLimits = [0.1 0.9];

inFlag = datenum(epiData.Date) <= fittedToDate;
tData = datenum(epiData.Date(inFlag))';
nCasesData = epiData.nCases(inFlag)';



% Setup earlyReject structure to disable early rejection for forward model
% (only used within ABC routine)
earlyReject.tData = [];
earlyReject.nCasesData = [];
earlyReject.threshold = inf;
earlyReject.rejectElimFlag = 1;

% Vector of dates for model simulations
t = par.date0 + (0:1:par.tEnd);

%  Sample from the posterior distribution, excluding particles whose value of Theta in the final time period is outside the specified percentiles
qt = quantile(Cvalues(:, end), ReffQuantLimits);        % find cut off points at specified quantiles
inFlag = Cvalues(:, end) >= qt(1) & Cvalues(:, end) <= qt(2);
Ctruncated = Cvalues(inFlag, :);                        % exclude Cvalues that are outside cut off points
CvaluesSample = Ctruncated( randsample( size(Ctruncated, 1), nReps, true), :);       % sample from truncated distribution

% Allocate arrays for storing model ouputs
nInfectedAll = zeros(nReps, par.tEnd+1);
nIsolAll = zeros(nReps, par.tEnd+1);
nIsolByAge = zeros(nReps, par.nAgeGroups);
nHospAll = zeros(nReps, par.tEnd+1);
nHospByAge = zeros(nReps, par.nAgeGroups);
nDiscAll = zeros(nReps, par.tEnd+1);
nDeathsAll = zeros(nReps, par.tEnd+1);
nDeathsByAge = zeros(nReps, par.nAgeGroups);
casesAll_0dose = zeros(nReps, par.tEnd+1);
casesAll_1dose = zeros(nReps, par.tEnd+1);
casesAll_2dose = zeros(nReps, par.tEnd+1);
ReffEmp = zeros(nReps, par.tEnd+1);
Rvt = zeros(nReps, par.tEnd+1);
Ct = zeros(nReps, par.tEnd+1);
TTIQeff_time = zeros(nReps, par.tEnd+1);
distRecent = zeros(nReps, 1);

% An array of parameter structures, one for each model realisation,
% initially all identical
parInd(1:nReps) = par;

% indices of the time vector to which trajectory ing will be applied
filterIndices =  (tData > fittedToDate-filterPeriod & tData <= fittedToDate);
tFilter = tData(filterIndices);
dataToFilter = nCasesData(filterIndices);

parfor iRep = 1:nReps
   parInd(iRep).Ct = makeCtVector( CvaluesSample(iRep, :), par);     %  modified parameter structure for the ith rep using the ith sample from the weighted posterior
   [cases, ~, Ct(iRep, :), ReffEmp(iRep, :), Rvt(iRep, :)] = runSimLeaky(parInd(iRep), earlyReject);     % run simulation model
   [nInf, nIsol, cases_0dose, cases_1dose, cases_2dose, nHosp, nDisc, ~, ~, nDeaths, TTIQeff_time(iRep, :) ] = postProcess(cases, parInd(iRep));      % post processing to extract model outputs over time
   nInfectedAll(iRep, :) = sum(nInf.');     % number infected by time summed over all age groups
   nIsolAll(iRep, :) = sum(nIsol.');        % number of confirmed cases by time summed over all age groups
   nIsolByAge(iRep, :) = sum(nIsol);        % number of confirmed cases by age group over all time       
   nHospAll(iRep, :) = sum(nHosp.');        % number of hospital admissions by time summed over all age groups
   nHospByAge(iRep, :) = sum(nHosp);        % number of confirmed cases by age group over all time       
   nDiscAll(iRep, :) = sum(nDisc.');        % number of hospital discharges by time summed over all age groups
   nDeathsAll(iRep, :) = sum(nDeaths.');    % number of deaths by time summed over all age groups
   nDeathsByAge(iRep, :) = sum(nDeaths);        % number of confirmed cases by age group over all time 
   casesAll_0dose(iRep, :) = sum(cases_0dose.'); 
   casesAll_1dose(iRep, :) = sum(cases_1dose.'); 
   casesAll_2dose(iRep, :) = sum(cases_2dose.');   
   distRecent(iRep) = calcError(tFilter, dataToFilter, t, nIsolAll(iRep, :));     % distance metric relative to the case data over the specified filtering period
end
keepFlag = distRecent < quantile(distRecent, filterFactor);     % keep only the realisations with the best fit to case data over the specified filtering period
CvaluesSample = CvaluesSample(keepFlag, :);
parInd = parInd(keepFlag);
ReffEmp = ReffEmp(keepFlag, :);
Rvt = Rvt(keepFlag, :);
Ct = Ct(keepFlag, :);
TTIQeff_time = TTIQeff_time(keepFlag, :);
nInfectedAll = nInfectedAll(keepFlag, :);
nIsolAll = nIsolAll(keepFlag, :);
nIsolByAge = nIsolByAge(keepFlag, :);
nHospAll = nHospAll(keepFlag, :);
nDiscAll = nDiscAll(keepFlag, :);
nDeathsAll = nDeathsAll(keepFlag, :);

nBedsAll = cumsum(nHospAll-nDiscAll, 2);        % hopsital beds occupied is cumulative admissions minus cumulative discharges

% Save results to file 
fOut = sprintf('results/results_fit_to_%s.mat', datestr(fittedToDate));
save(fOut);












