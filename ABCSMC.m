function [Theta, modelOutput, distances, weights] = ABCSMC(priorRnd, priorPDF, perturbRnd, perturbPDF, evalDist, modelOutSize, thresholds, nParticles, maxAcceptRate)

% USAGE: [Theta, modelOutput, distances, weights] = ABCSMC(priorRnd, priorPDF, perturbRnd, perturbPDF, evalDist, modelOutSize, thresholds, nParticles, maxAcceptRate)
% 
% NB: To run as simple ABC rejection, set threshold as scalar, priorPDF,
% perturnRnd and perturbPDF are not needed
%
% INPUTS: 
%    priorRnd -
%        function that when called as z = priorRnd() returns a 1 x nPars 
%        vector that is random samples from the joint prior distribution
%    priorPDF - 
%        function that when called as p = priorPDF(Theta), where Theta is
%        a 1 x nPars vector, returns the PDF of the prior distribution
%        evaluated at Theta
%    perturbRnd -
%        function that when called as y = perturbRnd(m) returns a m x nPars 
%        vector each of whose rows is a random sample from the perturbation 
%        distribution
%    perturbPDF -
%        function that when called as p = perturbPDF(Delta) where Delta is
%        a m x nPars vector returns a m x 1 column vector whose jth element
%        is the PDF of the perturbation distribution evaluated at the jth 
%        row of Delta = Theta1-Theta0
%    evalDist -
%        function that when called as [d, output] = evalDist(Theta, threshold)
%        where Theta is a 1 x nPars vectors simulates the forward model under
%        parameters  Theta and returns the distance metric relative to the 
%        data (d) and a vector of model output (output). Including a threshold 
%        in the call allows optional early termination of simulations in
%        which the threshold is exceed 
%    modelOutSize - 
%        scalar corresponding to the length of the vector that is output by
%        calls to evalDist
%    thresholds -
%        montonically decreasing row vector of acceptance thresholds for
%        each stage of the ABC-SMC algorithm
%    nParticles -
%        number of particles to retin at each stage of the algorithm
%    maxAcceptRate -
%        maximum acceptance rate at each stage of the ABCSMC algorithm. If
%        a higher proportion of particles than this meets the pre-specified
%        threshold, the alogirhtm will continue to draw candidates from the
%        proposal distribution and retain those with the smallest distance
%        metrics
%
% OUTPUTS
%    Theta - 
%        nParticles x nPar matrix whose rows represent the accepted particles
%    modelOutput - 
%        nParticles x modelOutSize matrix whose rows represent the 
%        corresponding model output that was returned for each particle
%    distances - 
%        nParticles x 1 vector of distances that were returned for each of
%        the accepted particles
%    weights -
%        nParticles x 1 vector of particle weights. The weighted 
%        distribution of particles contained in Theta is an ABC approximation
%        to the posterior distribution 

minSims = nParticles/maxAcceptRate;     % minimum number of simulations required to generate each population of particles

nThresholds = length(thresholds);       % number of populations required
nPars = length(priorRnd(1));            % number of parameters that are target for inference

Theta = nan(minSims, nPars);
distances = nan(minSims, 1);
modelOutput = nan(minSims, modelOutSize);
nAccepted = 0;
iSim = 0;
earlyRejThresh = thresholds(1);
while nAccepted < nParticles | iSim < minSims
    distProposed = zeros(nParticles, 1);
    outProposed = zeros(nParticles, modelOutSize);
    ThetaProposed = priorRnd(nParticles);  
    parfor iTrial = 1:nParticles
        [distProposed(iTrial), outProposed(iTrial, :)] = evalDist(ThetaProposed(iTrial, :), earlyRejThresh);
    end
    goodFlag = distProposed < earlyRejThresh;
    nGood = sum(goodFlag);
    Theta(nAccepted+1:nAccepted+nGood, :) = ThetaProposed(goodFlag, :);
    distances(nAccepted+1:nAccepted+nGood) = distProposed(goodFlag);
    modelOutput(nAccepted+1:nAccepted+nGood, :) = outProposed(goodFlag, :);
    nAccepted = nAccepted+nGood;
    if nAccepted > nParticles
       s = sort(distances(1:nAccepted));
       earlyRejThresh = s(nParticles);          % reduce threshold to the kth smallest found so far where k is the number of particles to be retained at the end
    end
    iSim = iSim+nParticles;
    fprintf('Stage %i/%i:  %i/%i accepted out of %i simulations (%.3f%%), threshold %.1f\n', 1, length(thresholds), nAccepted, nParticles, iSim, nAccepted/iSim*100, earlyRejThresh)
end
fprintf('\n');

% Retain best nParticles
[~, si] = sort(distances);
distances = distances(si(1:nParticles));
Theta = Theta(si(1:nParticles), :);
modelOutput = modelOutput(si(1:nParticles), :);
weights = 1/nParticles*ones(nParticles, 1);


for iThreshold = 2:nThresholds
    save(sprintf('interimResults%i.mat', iThreshold-1))
    
    ThetaPrevious = Theta;
    
    [~, popnSigma] = meanCovWeighted(ThetaPrevious, weights);       % calculate covariance matrix of previous population of particles
   
    Theta = nan(minSims, nPars);
    distances = nan(minSims, 1);
    modelOutput = nan(minSims, modelOutSize);
    nAccepted = 0;
    iSim = 0;
    earlyRejThresh = min(earlyRejThresh, thresholds(iThreshold));
    while nAccepted < nParticles | iSim < minSims
        distProposed = zeros(nParticles, 1);
        outProposed = zeros(nParticles, modelOutSize);
        proposedValid = zeros(nParticles, 1);
        ThetaProposed = nan(nParticles, nPars);
        while sum(proposedValid) < nParticles
            ind = (proposedValid == 0);
            ThetaProposed(ind, :) = ThetaPrevious( randsample(nParticles, sum(ind), true, weights), :) + perturbRnd(popnSigma, sum(ind)); 
            proposedValid(ind) = priorPDF(ThetaProposed(ind, :)) > 0;
        end        
        parfor iTrial = 1:nParticles
            [distProposed(iTrial), outProposed(iTrial, :)] = evalDist(ThetaProposed(iTrial, :), earlyRejThresh);
        end
        goodFlag = distProposed < earlyRejThresh;
        nGood = sum(goodFlag);
        Theta(nAccepted+1:nAccepted+nGood, :) = ThetaProposed(goodFlag, :);
        distances(nAccepted+1:nAccepted+nGood) = distProposed(goodFlag);
        modelOutput(nAccepted+1:nAccepted+nGood, :) = outProposed(goodFlag, :);
        nAccepted = nAccepted+nGood;
        if nAccepted > nParticles
           s = sort(distances(1:nAccepted));
           earlyRejThresh = s(nParticles);   % reduce threshold to the kth smallest found so far where k is the number of particles to be retained at the end
        end
        iSim = iSim+nParticles;
        fprintf('Stage %i/%i:  %i/%i accepted out of %i simulations (%.3f%%), threshold %.1f\n', iThreshold, nThresholds, nAccepted, nParticles, iSim, nAccepted/iSim*100, earlyRejThresh)
    end
    
    % Retain best nParticles
    [~, si] = sort(distances);
    distances = distances(si(1:nParticles));
    Theta = Theta(si(1:nParticles), :);
    modelOutput = modelOutput(si(1:nParticles), :);     
    
    % calculate weights
    Delta = repmat(Theta, nParticles, 1) - repelem(ThetaPrevious, nParticles, 1);
    ppdfMat = reshape( perturbPDF(Delta, popnSigma), nParticles, nParticles);
    y = priorPDF(Theta)./sum(weights.' .* ppdfMat, 2);
    weights = y/sum(y);
    ess = 1/sum(weights.^2)
end


