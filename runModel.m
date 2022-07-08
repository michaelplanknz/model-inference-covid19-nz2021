function [dist, nCasesSim] = runModel(Theta, tData, nCasesData, par, threshold)

par.Ct = makeCtVector(Theta, par);

earlyReject.tData = tData;
earlyReject.nCasesData = nCasesData;
earlyReject.threshold = threshold;

[cases, dist] = runSimLeaky(par, earlyReject);

tExt = par.date0 + (0:1:par.tEnd+1);
nCasesSim = histcounts(cases.tIsol, tExt);


