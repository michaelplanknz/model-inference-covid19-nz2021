clear all
close all

nDraws = 1000;
nToPlot = 20;
date0 = datenum('10AUG2021');
date1 = datenum('09JAN2022');



nPars = floor((date1-(date0+1)-14)/7);  

ta = date0+5 + 7*(0:nPars-1);

priorMean = [1.0, 0.5*ones(1, nPars-1)];           
priorSD =   [0.2, 0.5*ones(1, nPars-1)];           
corrLength = 2;                   

[ln_mu, ln_sigma] = logNormPars(priorMean, priorSD);                    % Get mu and sigma for a (series of univariate) lognormal distributions with the mean and SD specified above
priorSigma = zeros(nPars);
priorSigma(1, 1) = ln_sigma(1)^2;                                         % pre lockdown C(t) is independent of other values
priorSigma(2:end, 2:end) = getMvnSigma(ln_sigma(2:end), corrLength);      % Get a covariance matrix for a multivariate normal with the specified SD for each coordiate and correlation lengthscale corrLength

% Define prior and perturbation distributions

% when applying this you are working with parameters Theta in logged space
% - need to exponentiate in model before creating C(t) vector
priorRnd = @(nReps)(mvnrnd(ln_mu, priorSigma, nReps));  

Theta = priorRnd(nDraws);


qt = [0.1 0.25 0.75 0.9];      % quantiles to plot as uncertainty bands
clr = [0 0 1];
ls = '-';
mkr = 'none';


figure;
errorShade(ta, median(exp(Theta)), quantile(exp(Theta), qt),  clr, ls, mkr); 
hold on
plot(ta, exp(Theta(1:nToPlot, :)), 'Color', [0.4 0.4 0.4])
datetick('x', 'dd-mmm', 'keeplimits')
ylabel('control function C(t)')
xlim(ta([1 end]))


