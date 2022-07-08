function [cases, dist, Ct, ReffEmp, Rvt] = runSimLeaky(par, earlyReject)

% transmission process will stop if this many people have been infected
maxCases = 5100000;

% Set up time array
t = par.date0 + (0:1:par.tEnd);
nSteps = length(t)-1;

% Set the relative transmission rate due to alert level settings, as a function of time
% (this is the default setting which can be modified during the simulation by dynamic control rules if required but not implemented in this model)
Ct = par.Ct;

% calculate the area under the curve of the generation time distribution in 1-day time
% steps - this gives the relative amount of transmitting each person does
% on each day of their infection
tArr = 0:1:par.maxInfectTime+1;
C = wblcdf(tArr, par.genA, par.genB);       % CDF of the generation time distribution
auc = diff(C).';                            % area under the generation time distribution in one day increments
auc = auc/sum(auc);                         % renormalise auc so that it sums to 1 


% initialise variables for the case table
caseID = (1:maxCases).';
parentID = nan(maxCases, 1);
gen = nan(maxCases, 1);
nOff = zeros(maxCases, 1);
ageGroup = nan(maxCases, 1);
Rimult = genRi(maxCases, 1, par);            % pre generate each individual's values of Yi
subclinFlag = nan(maxCases, 1);
vaccDoses = zeros(maxCases, 1);
tInfect = nan(maxCases, 1);
tOnset = genOnsetDelay(maxCases, 1, par);   % pre-generate each individuals incubation period
tIsol = nan(maxCases, 1);                   % testing and isolation times assumed simulteneous
tQuar = nan(maxCases, 1);                    % testing and isolation times assumed simulteneous
tHosp = nan(maxCases, 1);
tDisc = nan(maxCases, 1); 
icuFlag = zeros(maxCases, 1); 
diedFlag = zeros(maxCases, 1); 

% initialise case table 
cases = table(caseID, parentID, gen, nOff, ageGroup, Rimult, subclinFlag, vaccDoses, tInfect, tOnset, tQuar, tIsol, tHosp, tDisc, icuFlag, diedFlag);

% define number of seed cases and their properties
nSeedCases = poissrnd(par.meanSeedsPerDay*par.seedPeriod) + 0;
cases.parentID(1:nSeedCases) = 0;
cases.gen(1:nSeedCases) = 1;
cases.ageGroup(1:nSeedCases) = randsample( length(par.popDist), nSeedCases, true, par.popDist);    % ages of seed cases distributed the same as the overall popn
cases.Rimult(1:nSeedCases) = par.relInfSeedCases * cases.Rimult(1:nSeedCases);                     % seed cases can have an infectiousness multiplier, e.g. to represent home isolation of border cases
cases.vaccDoses(1:nSeedCases) = par.genSeedVaxStatus(cases, nSeedCases, par);
cases.subclinFlag(1:nSeedCases) = rand(nSeedCases, 1) < par.pSub(cases.ageGroup(1:nSeedCases));
cases.tInfect(1:nSeedCases) = t(1) + floor(par.seedPeriod*rand(nSeedCases, 1));
cases.tOnset(1:nSeedCases) = cases.tOnset(1:nSeedCases) + cases.tInfect(1:nSeedCases);
pIsol = par.pTestClin*(cases.subclinFlag(1:nSeedCases) == 0 ) + par.pTestSub*(cases.subclinFlag(1:nSeedCases) == 1 );
isolFlag = rand(nSeedCases, 1) < pIsol;
nIsol = sum(isolFlag);
tIsol = cases.tOnset(isolFlag) + genIsolDelay(nIsol, 1, par);
ind = tIsol < t(1)+par.minDetectTime;
tIsol(ind) = t(1)+par.minDetectTime + rand(sum(ind), 1)*par.followUpTime;
cases.tIsol(isolFlag) = tIsol;

% Fractions of each age group that is in the susceptible class and has had
% 0, 1 or 2 doses of vaccines. These will initially sum to 1 but will
% become depleted as people get infected
susFrac0 = zeros(nSteps, par.nAgeGroups);
susFrac1 = zeros(nSteps, par.nAgeGroups);
susFrac2 = zeros(nSteps, par.nAgeGroups);
susFrac0(1, :) = 1-par.cov1(1, :)-par.cov2(1, :);
susFrac1(1, :) = par.cov1(1, :);
susFrac2(1, :) = par.cov2(1, :);

% Fractions moving between dose compartments each time step
q01 = min(1, max(0, (par.cov1(2:end, :)+par.cov2(2:end, :)-par.cov1(1:end-1, :)-par.cov2(1:end-1, :))./(1-par.cov1(1:end-1, :)-par.cov2(1:end-1, :))));
q01 = [q01; zeros(1, par.nAgeGroups)];
q12 = min(1, max(0, (par.cov2(2:end, :)-par.cov2(1:end-1, :))./par.cov1(1:end-1, :)));
q12 = [q12; zeros(1, par.nAgeGroups)];

% Relative transmisison rates as a result of vaccination and
% quarantine/isolation
relTransVacc = [1; 1-par.VEt1; 1-par.VEt2];
relTransIsol = [1; par.cQuar; par.cIsol];


nCases0 = histcounts(cases.ageGroup(1:nSeedCases), 1:par.nAgeGroups+1);
nCases = nCases0;
ReffEmp = zeros(size(t));
Rvt = zeros(size(t));
dailyAvg = 0;
nActive = 0;
nFuture = nSeedCases;
iStep = 1;
dist = 0;

while iStep < nSteps & sum(nCases) < maxCases & dist < earlyReject.threshold 
   
    iStep = iStep+1;
    
    % number of seed cases whose infeciotn time is still in the future:
    nFuture = sum( t(iStep) <= cases.tInfect );
    
    % Get IDs of active cases at current time step:
    activeID = cases.caseID(t(iStep) > cases.tInfect & t(iStep) <= cases.tInfect+par.maxInfectTime);
    nActive = length(activeID);

    if nActive > 0
        % Simulate new secondary cases on day i from each currently active
        % case:
        
        % Area under the curve of the transmission rate for each active case at current time step
        auci = auc(t(iStep)-cases.tInfect(activeID));       
        
        % isolation status of active cases (0 for nothing, 1 for
        % quarantined, 2 for isolated)
        isolStatus = (t(iStep) > cases.tQuar(activeID) & ~(t(iStep) > cases.tIsol(activeID))) + 2*(t(iStep) > cases.tIsol(activeID));
        
        % Calculate the expected number of offspring in each age group during the current
        % time step from each active case, assuming a fully susceptible population.
        % This defines a matrix expOff whose i,j element is the expected number of
        % offspring from parent case i in age group j
        % This is the product of the following factors for each case:
        % - relative transmission due to current alert level
        % - relative transmission due to vaccination status of active case i
        % - relative reproduction number (typically gamma distributed) of case i
        % - relative transmission rate for clinical/subclinical status of case i
        % - relative transmission due to isolation/quarantine status of case i
        % - time-dependent transmission rate for the parent case i, quantified by auci
        % - jth column of the NGM for an unvaccinated clinical individual in the age group of the parent case
        expOff = Ct(iStep) * relTransVacc(1+cases.vaccDoses(activeID)) .* relTransIsol(1+isolStatus) .* cases.Rimult(activeID) .* (1 - (1-par.cSub)*cases.subclinFlag(activeID)) .* auci .* (par.NGMclin(:, cases.ageGroup(activeID))).';
        
        % split expOff into expected offspring who have had 0, 1 or 2 doses
        % In doing this, some putative infections are prevented by the
        % vaccine (and others by immunity from prior infection is susFrac0,
        % susFrac1 and susFrac2 sum to less than 1).
        expOff0 = expOff .* susFrac0(iStep-1, :);
        expOff1 = (1-par.VEi1) * expOff .* susFrac1(iStep-1, :);
        expOff2 = (1-par.VEi2) * expOff .* susFrac2(iStep-1, :);
        ReffEmp(iStep) = sum(sum(expOff0+expOff1+expOff2))/sum(auci);
        
        % Generate nOff actual number of offspring from parent case i in age group j, thinned due to immunity in each age group as a result of infection prevention in vaccinated individuals and prior infection prevention in vaccinated offspring
        nOff0 = poissrnd(expOff0);   
        nOff1 = poissrnd(expOff1); 
        nOff2 = poissrnd(expOff2); 
        
        % Total number of offspring summed across all parent cases and all
        % age groups:
        nOffTot = sum(sum(nOff0+nOff1+nOff2));

        if nOffTot > 0
            secIDs = (sum(nCases)+1:sum(nCases)+nOffTot).';       % IDs for today's newly infected cases

            % assign age groups, parent ID, clinical status and vaccination status for new
            % cases based on nOff matrices
            ageGroupList = repmat(1:par.nAgeGroups, nActive, 1);
            parentNumber = repmat(1:nActive, 1, par.nAgeGroups);
            nx = par.nAgeGroups*nActive;
            % Make an array (propList) whose coluumns are: (1) age group of
            % offspring, (2) ID of parent case, (3) number of doses of offspring
            propList = [repmat([ageGroupList(:), parentNumber(:)], 3, 1), [zeros(nx, 1); ones(nx, 1); 2*ones(nx, 1)] ];
            % The frequency (number of cases with each age gorup, parent ID
            % and dose combination) of each row of propList is given by
            % the elements in the nOff matrices. Generate a sample X whose
            % columns are the values of these three properties for each
            % secondary case:
            X = repelem( propList, [nOff0(:); nOff1(:); nOff2(:)], 1);
            cases.ageGroup(secIDs) = X(:, 1);
            cases.parentID(secIDs) = activeID( X(:, 2) );
            cases.vaccDoses(secIDs) = X(:, 3);
            v0Flag = cases.vaccDoses(secIDs) == 0;
            v1Flag = cases.vaccDoses(secIDs) == 1;
            v2Flag = cases.vaccDoses(secIDs) == 2;
            cases.subclinFlag( secIDs(v0Flag) ) = rand(sum(v0Flag), 1) < par.pSub(cases.ageGroup(secIDs(v0Flag)));
            cases.subclinFlag( secIDs(v1Flag) ) = rand(sum(v1Flag), 1) < par.VEs1 + (1-par.VEs1)*par.pSub(cases.ageGroup(secIDs(v1Flag)));
            cases.subclinFlag( secIDs(v2Flag) ) = rand(sum(v2Flag), 1) < par.VEs2 + (1-par.VEs2)*par.pSub(cases.ageGroup(secIDs(v2Flag)));
            
            cases.gen(secIDs) = cases.gen(cases.parentID(secIDs))+1;    % Generation of each new cases is generation of parent + 1
            cases.tInfect(secIDs) = t(iStep);                           % Infection time for each new cases is today
            cases.tOnset(secIDs) = t(iStep) + cases.tOnset(secIDs);     % For efficiency infection to onset delay is pre-stored in cases.tOnset
            cases.nOff(activeID) = cases.nOff(activeID)+sum(nOff0+nOff1+nOff2, 2);      

            % simulate case testing and isolation effects for new cases
            pIsol = par.pTestClin*(cases.subclinFlag(secIDs) == 0 ) + par.pTestSub*(cases.subclinFlag(secIDs) == 1 ) ;
            isolFlag = rand(nOffTot, 1) < pIsol;
            nIsol = sum(isolFlag);
            tIsol = cases.tOnset(secIDs(isolFlag)) + genIsolDelay(nIsol, 1, par);
            ind = tIsol < t(1)+par.minDetectTime;
            tIsol(ind) = t(1)+par.minDetectTime + rand(sum(ind), 1)*par.followUpTime;
            cases.tIsol(secIDs(isolFlag)) = tIsol;
            
            % simulate contact tracing effects for new cases
            pTrace = par.pTrace * (~(dailyAvg > par.traceCapacity)) *(~isnan(cases.tIsol(cases.parentID(secIDs))));
            traceFlag = rand(nOffTot, 1) < pTrace;
            nTrace = sum(traceFlag);
            cases.tQuar(secIDs(traceFlag)) = cases.tIsol(cases.parentID(secIDs(traceFlag))) + genTraceDelay(nTrace, 1, par);
            % optional: individuals traced prior to onset go into full isolation (as opposed to quarantine) on symptom onset:
            % this is applied to all individuals including subclinical on
            % assumption that asymtomatic contacts get tested. This allows
            % offspring of subclinicals to be traced
            cases.tIsol(secIDs(traceFlag)) = min(cases.tIsol(secIDs(traceFlag)), max(cases.tQuar(secIDs(traceFlag)), cases.tOnset(secIDs(traceFlag))));

            % simulate hospitalisation outcomes for new cases (note
            % hospitalisations are done here because hospitalisation leads to
            % detection if case has not already tested; deaths are done after 
            % the simulation is complete as they don't affect epidemic dynamics)
            clinFlag = cases.subclinFlag(secIDs) == 0;
            pHosp = ((cases.vaccDoses(secIDs) == 0) + (1-par.VEd1)/(1-par.VEs1)*(cases.vaccDoses(secIDs) == 1) + (1-par.VEd2)/(1-par.VEs2)*(cases.vaccDoses(secIDs) == 2)) .* par.IHR(cases.ageGroup(secIDs))./par.IDR(cases.ageGroup(secIDs));
            hospFlag = clinFlag & (rand(nOffTot, 1) < pHosp);
            cases.tHosp(secIDs(hospFlag)) = cases.tOnset(secIDs(hospFlag)) + genHospDelay(sum(hospFlag), 1, par);
            cases.tDisc(secIDs(hospFlag)) = cases.tHosp(secIDs(hospFlag)) + genHospLOS(sum(hospFlag), 1, par);
            % cases are detected and isolated once hospitalised
            cases.tIsol(secIDs(hospFlag)) = min(cases.tIsol(secIDs(hospFlag)), cases.tHosp(secIDs(hospFlag)));
           
            % Update cumulative infections to date (in each age group)
            nCases = nCases+sum(nOff0+nOff1+nOff2, 1);        
        end
    else
        nOff0 = zeros(1, par.nAgeGroups);
        nOff1 = zeros(1, par.nAgeGroups);
        nOff2 = zeros(1, par.nAgeGroups);
    end   

    % Update susceptible fractions according to vaccinations given
    % and new infections this time step
    susFrac2(iStep, :) = max(0, (susFrac2(iStep-1, :)  - sum(nOff2, 1)./(par.popCount')) + (susFrac1(iStep-1, :)- sum(nOff1, 1)./(par.popCount')).*q12(iStep-1, :) );
    susFrac1(iStep, :) = max(0, (susFrac1(iStep-1, :)- sum(nOff1, 1)./(par.popCount')).*(1-q12(iStep-1, :)) + (susFrac0(iStep-1, :)- sum(nOff0, 1)./(par.popCount')).*q01(iStep-1, :) ) ;
    susFrac0(iStep, :) = max(0, (susFrac0(iStep-1, :)- sum(nOff0, 1)./(par.popCount')).*(1-q01(iStep-1, :))) ;
   
    infBlock = susFrac0(iStep-1, :) + (1-par.VEi1)*susFrac1(iStep-1, :) + (1-par.VEi2)*susFrac2(iStep-1, :);    % row
    transBlockNum0 =                               (par.IDR + par.cSub*(1-par.IDR))                           .* susFrac0(iStep-1, :)';
    transBlockNum1 = (1-par.VEi1) * (1-par.VEt1) * ((1-par.VEs1)*par.IDR + par.cSub*(1-(1-par.VEs1)*par.IDR)) .* susFrac1(iStep-1, :)';
    transBlockNum2 = (1-par.VEi2) * (1-par.VEt2) * ((1-par.VEs2)*par.IDR + par.cSub*(1-(1-par.VEs2)*par.IDR)) .* susFrac2(iStep-1, :)';  
    transBlock = (transBlockNum0+transBlockNum1+transBlockNum2)./infBlock';                                     % column
    NGMt = par.NGMclin .* infBlock .* transBlock;
    Rvt(iStep) = eigs(NGMt, 1);
    
    
    % Update 7-day average of new daily cases
    dailyAvg = sum( t(iStep) > min(cases.tIsol, cases.tQuar) & t(iStep) <= min(cases.tIsol, cases.tQuar)+7 )/7;
    
    % calculate distance metric to allow early rejection of simulation
    % during fitting process
    nIsolTemp = histcounts(cases.tIsol, [t(1):t(iStep)+1] );
    dist = calcError(earlyReject.tData, earlyReject.nCasesData, t(1):t(iStep), nIsolTemp);
end

% If simulation stopped early, assign a sufficiently large value to the
% distance metric to ensure the simulation is rejected by the fitting
% algorithm
if iStep < nSteps
   dist = 2*earlyReject.threshold; 
end

totCases = sum(nCases);
cases = cases(1:totCases, :);
cases.icuFlag = ~isnan(cases.tHosp) & rand(totCases, 1) < par.pICU(cases.ageGroup)  ;  
cases.diedFlag = ~isnan(cases.tHosp) & rand(totCases, 1) < par.IFR(cases.ageGroup)./par.IHR(cases.ageGroup);  






