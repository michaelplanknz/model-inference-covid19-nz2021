function x = myPerturbRnd( priorSigma, popnSigma, nReps)

% Generate deviates from the perturbation distribution

pertSigma = getPertSigma(priorSigma, popnSigma);

[nPars, ~] = size(priorSigma);
x = mvnrnd(zeros(1, nPars), pertSigma, nReps);


