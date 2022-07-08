function par = getPar(vaxDataDate)


%------------- Seed and control Parameters --------------
par.R0 = 6;  
par.tEnd = 174;  
par.date0 = datenum('10AUG2021');

par.Ct = ones(1, par.tEnd+1);

par.minDetectTime = 7;    % time of outbreak detection 
par.followUpTime = 7;    % cases with an isolation time prior to detection are distributed over this time period post-detection
par.seedPeriod = 8;
par.meanSeedsPerDay = 200/par.seedPeriod;       
par.relInfSeedCases = 1;
par.genSeedVaxStatus = @genSeedVaxStatus_copyPopn;


%------------- Branching Process Parameters --------------

par.cSub = 0.5; % Relative infectiousness of subclinicals
par.ssk = 0.5; % Overdispersion/superspreading parameter k
par.maxInfectTime = 21; % Maximum infectious perod 
par.genA = 5.665; par.genB = 2.826; % Generation time distribution parameters
par.incA = 5.8; par.incB = 0.95; % Exposure to Onset distribution parameters

par.isolA = 1; par.isolB = 4; % Onset to isolation distribution parameters
par.traceA = 3; par.traceB = 1; % Parent isolation to quarantine distribution parameters
par.hospA = 1; par.hospB = 5; % Onset to hospitalisation distribution parameters
par.losA = 1; par.losB = 8; % Hospital LOS distribution parameters

par.pTestClin = 0.45; par.pTestSub = 0; % Probability of detecting symptomatic / subclinical case  (45% for both)
par.pTrace = 0.7; % Probability of detecting a case by contact tracing
par.traceCapacity = inf;        % contact tracing capacity (average daily cases over last 7 days)

par.cIsol = 0;
par.cQuar = 0.5;


%------------- Vaccine Effectiveness Parameters --------------

par.VEi1 = 0.55; % Default vaccine effectiveness against infection 
par.VEs1 = 0;    % 0
par.VEt1 = 0;    % 0  Default vaccine effectiveness against transmission given infection
par.VEd1 = 0.6; % Default vaccine effectiveness against severe disease (and death) given infection - change to 0.6 to roughly agree with PHE data

par.VEi2 = 0.7; % Default vaccine effectiveness against infection
par.VEs2 = 0.5;   % 0    [0.5] Vaccine effctiveness against symptomatic disease given infection.
par.VEt2 = 0.375; % 0.5  [0.375]  Default vaccine effectiveness against transmission given infection. 0.375 if including VEs = 0.5 for reduction in symptomatic breakthroughs
par.VEd2 = 0.8; % Default vaccine effectiveness against severe disease (and death) given infection




%------------- Disease Rate Data --------------
par.nAgeGroups = 16;
par.IDR = [0.5440, 0.5550, 0.5770, 0.5985, 0.6195, 0.6395, 0.6585, 0.6770, 0.6950, 0.7117, 0.7272, 0.7418, 0.7552, 0.7680, 0.7800, 0.8008]';    % probability of symptomatic infection by age
par.ui = [0.4000, 0.3950, 0.3850, 0.4825, 0.6875, 0.8075, 0.8425, 0.8450, 0.8150, 0.8050, 0.8150, 0.8350, 0.8650, 0.8450, 0.7750, 0.7400 ];    % relative suscewptibility by age
par.pSub = 1-par.IDR;

[par.IHR, par.pICU, par.IFR] = getHerreraRates();



%------------- Specify Population Structure --------------
fNamePop = 'data/nzpopdist.xlsx';
fprintf('   Loading population distribution:    %s\n', fNamePop)
popSizeData = readmatrix(fNamePop); % Load NZ population structure from data folder

par.popCount = zeros(par.nAgeGroups, 1); % Create popDist vector
par.popCount(1:par.nAgeGroups-1) = popSizeData(1:par.nAgeGroups-1, 2); % Fill entries with population distribution
par.popCount(par.nAgeGroups) = sum(popSizeData(par.nAgeGroups:end, 2)); % Aggregate 75+ age-groups
par.totalPopSize = sum(par.popCount);
par.popDist = par.popCount/sum(par.popCount); 

par.age = (2.5:5:77.5)'; % Define final age groups (matching contact matrix)


%------------- Get next generation matrix ---------
[par.NGM, par.NGMclin] = getNGM(par);



%------------- Define vaccine coverage either static or time-dependent ---------
par.vaccImmDelay = 14;     % delay from vaccination to immunity
convertPeriod = 56;                        % dose conversions will be applied over this time period (days)


% Creates matrices for 1st and 2nd dose coverage by time and age group, starting
% on date0
% If the simulation runs longer than the number of rows in these matrices, the
% coverage on the last defined time date will be used for the remainder of 
% the simulation
[dates, vacc1, vacc2] = getVaccineData(vaxDataDate);    % read time-dependent vaccination data
[vacc1, vacc2] = convertFirstDoses(datenum(dates), vacc1, vacc2, vaxDataDate+1, convertPeriod, par);
dates = dates+par.vaccImmDelay;     % shift vaccination dates to allow for delay
iDate = find(datenum(dates) == par.date0);
nPad = max(0, par.tEnd - (length(dates) - iDate));      % pad vaccine data if necessary by assuming coverage remains fixed at final data point
dates = [dates(iDate:end), dates(end)+1:dates(end)+nPad];
par.cov1 = [vacc1(iDate:end, :); repmat(vacc1(end, :), nPad, 1)];
par.cov2 = [vacc2(iDate:end, :); repmat(vacc2(end, :), nPad, 1)];
[nVaxDates, ~] = size(par.cov1);


cov1_12plus = sum(par.cov1(:, 3:end).*par.popCount(3:end).', 2)/( 3/5*par.popCount(3) + sum(par.popCount(4:end)) );
cov2_12plus = sum(par.cov2(:, 3:end).*par.popCount(3:end).', 2)/( 3/5*par.popCount(3) + sum(par.popCount(4:end)) );


% get dominant eigenvalue of vaccinated NGM through time to calculater
% effect of vaccination on reproduciton number
Rvi = zeros(1, nVaxDates);
for iDate = 1:nVaxDates
    Rvi(iDate) = calcVaxedR(1-par.cov1(iDate, :)-par.cov2(iDate, :), par.cov1(iDate, :), par.cov2(iDate, :), par);
end


