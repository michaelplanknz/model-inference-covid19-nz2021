function [Sigma] = getMvnSigma(priorSD, corrLength)

n = length(priorSD);
x = linspace(1, n, n)'; % col vector of centres of locations for Reff values

C = exp(-(x-x').^2/(2*corrLength^2));
D = diag(priorSD);
Sigma = D*C*D;

