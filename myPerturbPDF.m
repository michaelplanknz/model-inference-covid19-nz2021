function f = myPerturbPDF(x, priorSigma, popnSigma)

% Evaluate PDF of perturbation distribution

pertSigma = getPertSigma(priorSigma, popnSigma);

[nPars, ~] = size(priorSigma);
f = mvnpdf(x, zeros(1, nPars), pertSigma);



