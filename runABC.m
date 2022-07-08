function runABC(epiData, fittedToDate, par)


nParticles = 1000;       % 1000
maxAcceptRate = 0.01;     % 0.01
thresholds = 6500*0.8.^[0, 1, 2];        % 450 for full to 30 Oct, 700 to 3 Nov

nPars = floor((fittedToDate-(par.date0+1)-14)/7);       % Leave at least 2 weeks between the end of the last fitted period for C(t) and the last day of data 
priorMean = [1.0, 0.5*ones(1, nPars-1)];            % was 0.75
priorSD =   [0.2, 0.5*ones(1, nPars-1)];            % was 0.75
corrLength = 2;                     % WAS 2 until March 2022 - correlation length in terms of number of time periods for C(t): 1/l how fast correlation drops off

[ln_mu, ln_sigma] = logNormPars(priorMean, priorSD);                    % Get mu and sigma for a (series of univariate) lognormal distributions with the mean and SD specified above
priorSigma = zeros(nPars);
priorSigma(1, 1) = ln_sigma(1)^2;                                         % pre lockdown C(t) is independent of other values
priorSigma(2:end, 2:end) = getMvnSigma(ln_sigma(2:end), corrLength);      % Get a covariance matrix for a multivariate normal with the specified SD for each coordiate and correlation lengthscale corrLength

% Define prior and perturbation distributions

% when applying this you are working with parameters Theta in logged space
% - need to exponentiate in model before creating C(t) vector
priorRnd = @(nReps)(mvnrnd(ln_mu, priorSigma, nReps));  
priorPDF = @(Theta)(mvnpdf(Theta, ln_mu, priorSigma));  

perturbRnd = @(popnSigma, nReps)myPerturbRnd(priorSigma, popnSigma, nReps);
perturbPDF = @(Delta, popnSigma)myPerturbPDF(Delta, priorSigma, popnSigma);


inFlag = datenum(epiData.Date) <= fittedToDate;
tData = datenum(epiData.Date(inFlag))';
nCasesData = epiData.nCases(inFlag)';

% Run ABCSMC algorithm
parABC = par;
parABC.tEnd = tData(end)-par.date0;            % modify model runtime so that runs up to the specified fittedToDate
modelOutSize = parABC.tEnd+1;

% Define function that returns the distance metric between model output and
% data:
% Note: exponentiate Theta before passing to simulation model:
evalDist = @(lnTheta, threshold)(runModel(exp(lnTheta), tData, nCasesData, parABC, threshold));

[lnTheta, modelOutput, distances, weights] = ABCSMC(priorRnd, priorPDF, perturbRnd, perturbPDF, evalDist, modelOutSize, thresholds, nParticles, maxAcceptRate);
%%
% Save results to file
save(sprintf('results/ABCSMC_cases_to_%s.mat', datestr(fittedToDate)), 'lnTheta' ,'modelOutput', 'distances', 'weights'  )



figure
plot(tData, nCasesData, 'o', parABC.date0:parABC.date0+parABC.tEnd, quantile(modelOutput, [0.025,0.25,0.5,0.75,0.975]))
datetick('x', 'dd-mmm', 'keeplimits')
ylabel('daily reported cases')

lnThetaSample = lnTheta(randsample(length(weights), 1e5, true, weights), :);


figure
qtl = quantile(exp(lnThetaSample), [0.25 0.5 0.75]);
dqtl = diff(qtl);
errorbar(1:nPars, qtl(2, :), dqtl(1, :), dqtl(2, :)) 



