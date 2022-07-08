function [nInfected, nIsol, cases_0dose, cases_1dose, cases_2dose, nHosp, nDisc, nICUIn, nICUOut, nDeaths, TTIQeff_time] = postProcess(cases, par)

% Set up time array
t = par.date0 + (0:1:par.tEnd);

% Only count Reff over cases who have completed their infectious period
%Reff = mean(cases.nOff(cases.tInfect < t(iStep) - par.maxInfectTime));

% calculate TTIQ effect overall and over time
wI = 1-wblcdf(floor(cases.tIsol)-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if isolation was perfect
wQ = 1-wblcdf(floor(min(cases.tQuar, cases.tIsol))-cases.tInfect, par.genA, par.genB); % amount of transmission prevented if quarantine was perfect    % wQ >= wI
wI(isnan(wI)) = 0;
wQ(isnan(wQ)) = 0;
ReffReduc = (1-par.cQuar)*(wQ-wI) + (1-par.cIsol)*(wI);
infAfterDetectFlag = cases.tInfect >= t(1)+par.minDetectTime;
indSub = infAfterDetectFlag & cases.subclinFlag == 0;
indClin = infAfterDetectFlag & cases.subclinFlag == 1;
TTIQeff = 1 - (sum(1-ReffReduc(indClin)) + par.cSub*sum(1-ReffReduc(indSub))) / (sum(indClin) + par.cSub*sum(indSub));
TTIQeff_time = zeros(size(t));
for ii = 1:length(t)
    indSub = infAfterDetectFlag & cases.subclinFlag == 0 & (cases.tInfect <= t(ii) & cases.tInfect > t(ii)-7);
    indClin = infAfterDetectFlag & cases.subclinFlag == 1 & (cases.tInfect <= t(ii) & cases.tInfect > t(ii)-7);
    TTIQeff_time(ii) = 1 - (sum(1-ReffReduc(indClin)) + par.cSub*sum(1-ReffReduc(indSub))) / (sum(indClin) + par.cSub*sum(indSub));
end
TTIQeff_time(isnan(TTIQeff_time)) = 0;

% Count number of cases (and other outcomes) in each age group on each day
tExt = [t t(end)+1];
nInfected = histcounts2(cases.tInfect, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nIsol = histcounts2(cases.tIsol, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
cases_0dose = histcounts2(cases.tIsol(cases.vaccDoses == 0), cases.ageGroup(cases.vaccDoses == 0), tExt, 1:par.nAgeGroups+1);
cases_1dose = histcounts2(cases.tIsol(cases.vaccDoses == 1), cases.ageGroup(cases.vaccDoses == 1), tExt, 1:par.nAgeGroups+1);
cases_2dose = histcounts2(cases.tIsol(cases.vaccDoses == 2), cases.ageGroup(cases.vaccDoses == 2), tExt, 1:par.nAgeGroups+1);
nHosp = histcounts2(cases.tHosp, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nDisc = histcounts2(cases.tDisc, cases.ageGroup, tExt, 1:par.nAgeGroups+1);
nICUIn = histcounts2(cases.tHosp(cases.icuFlag == 1), cases.ageGroup(cases.icuFlag == 1), tExt, 1:par.nAgeGroups+1);
nICUOut = histcounts2(cases.tDisc(cases.icuFlag == 1), cases.ageGroup(cases.icuFlag == 1), tExt, 1:par.nAgeGroups+1);
nDeaths = histcounts2(cases.tDisc(cases.diedFlag == 1), cases.ageGroup(cases.diedFlag == 1), tExt, 1:par.nAgeGroups+1);





