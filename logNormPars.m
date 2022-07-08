function [Mu, Sigma] = logNormPars(m, sd)

% Calculates the parameters of a log normal distribution with mean m and
% standard deviation sd, suitable for use in Matlab logn functions.
%
% USAGE: [Mu, Sigma] = logNormPars(m, sd)
%
% INPUTS: m - mean of the required distribution
%         sd - standard deviation of the required distribution
%
% OUTPUTS: Mu, Sigma - parameters defining the appropriate log normal distribution 

Mu = log( m.^2 ./ sqrt(m.^2 + sd.^2));
Sigma = sqrt( log( 1 + sd.^2./m.^2 ));





