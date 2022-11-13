function fitToData(epiData, fittedToDate, par, hyperPar)


nPars = floor((fittedToDate-(par.date0+1)-14)/7);       % Leave at least 2 weeks between the end of the last fitted period for C(t) and the last day of data 

% Specify mean, SD and correlation lengthscale of prior for control function
priorMean = [1.0, 0.5*ones(1, nPars-1)];           
priorSD =   [0.2, 0.5*ones(1, nPars-1)];           
corrLength = 2;                                        % correlation lengthscale in terms of number of time periods for C(t): 1/l how fast correlation drops off

[ln_mu, ln_sigma] = logNormPars(priorMean, priorSD);                    % Get mu and sigma for a (series of univariate) lognormal distributions with the mean and SD specified above
priorSigma = zeros(nPars);
priorSigma(1, 1) = ln_sigma(1)^2;                                         % pre lockdown C(t) is independent of other values
priorSigma(2:end, 2:end) = getMvnSigma(ln_sigma(2:end), corrLength);      % Get a covariance matrix for a multivariate normal with the specified SD for each coordiate and correlation lengthscale corrLength
L = cholcov(priorSigma);                                                % Calculate Cholevksy decomposition will will be used to transform standard normals to normals with covariance matrix priorSigma

% Define prior and perturbation distributions (standard multivariate
% normals)
priorRnd = @(nReps)(randn(nReps, nPars));                       % faster to use randn than mvnrnd for normal random deviates with identity covariance matrix
priorLogPDF = @(Theta)( sum(normLogPDF(Theta, 0, 1), 2) );      % log PDF for of a multivariate normal with identity covariance matrix evaluated at a given row vector is just the sum of the univariate log PDF evaluated at each element of the vector 
perturbRnd = @(Sigma, nReps)(mvnrnd( zeros(1, size(Sigma, 1)), Sigma, nReps));

inFlag = datenum(epiData.Date) <= fittedToDate;
tData = datenum(epiData.Date(inFlag))';
nCasesData = epiData.nCases(inFlag)';
par.tEnd = tData(end)-par.date0;            % modify model runtime so that runs up to the specified fittedToDate


% Define function that returns the distance metric between model output and
% data:
% Note: exp(ln_mu + Theta*L) transforms deviates from an uncorrelated normal to a multivariate normal with covariance matrix prioSigma and then exponentiates before passing to simulation model
evalDist = @(Theta, threshold)(runModel(exp(ln_mu + Theta*L), tData, nCasesData, par, threshold));


% Run ABCSMC algorithm
particles = doABCSMC(priorRnd, priorLogPDF, perturbRnd, evalDist, hyperPar, false);
Theta = particles.Theta;
Cvalues = exp(ln_mu + Theta*L);
modelOutput = particles.modelOutput;
distances = particles.dist;


% Save results to file
save(sprintf('results/ABCSMC_cases_to_%s.mat', datestr(fittedToDate)), 'Theta', 'Cvalues', 'modelOutput', 'distances')



figure
plot(tData, nCasesData, 'o', par.date0:par.date0+par.tEnd, quantile(modelOutput, [0.025,0.25,0.5,0.75,0.975]))
datetick('x', 'dd-mmm', 'keeplimits')
ylabel('daily reported cases')


figure
qtl = quantile(Cvalues, [0.25 0.5 0.75]);
dqtl = diff(qtl);
errorbar(1:nPars, qtl(2, :), dqtl(1, :), dqtl(2, :)) 



