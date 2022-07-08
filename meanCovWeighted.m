function [Mu, Sigma] = meanCovWeighted(X, weights)

% Function calculate the mean and covariance matrix of weighted sample of m
% deviates in R^n 
%
% USAGE: [Mu, Sigma] = meanCovWeighted(X, weights)
%
% INPUTS: X       - m x n matrix whose ith row is the ith observation in the sample 
%         weights - m x 1 vector of weights
%
% OUTPUTS: Mu - 1 x n vector that is the weighted mean of the sample
%          Sigma - n x n weighted covariance matrix


nSample = length(weights);

weights = weights/sum(weights);     % normalise so weights sum to 1

Mu = sum(weights.*X);
Sigma = nSample/(nSample-1) * (weights.*(X-Mu))'*(X-Mu);


