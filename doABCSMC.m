function particles = doABCSMC(priorRnd, priorLogPDF, perturbRnd, evalDist, hyperPar, loadRejectStepFlag)

% ABC-SMC algorithm based on Carr, Drovandi and Simpson, Interface 2021
%
% USAGE: particles = doABCSMC(priorRnd, priorLogPDF, perturbRnd, evalDist, hyperPar, loadRejectStepFlag)
%
% INPUTS:
%    priorRnd -
%        function that when called as z = priorRnd(n) returns a n x nPars 
%        vector that is random samples from the joint prior distribution
%    priorLogPDF - 
%        function that when called as p = priorLogPDF(Theta), where Theta is
%        a m x nPars vector, returns a m x 1 vector containing the log PDF
%        of the prior distribution evaluated at each row of Theta
%    perturbRnd -
%        function that when called as y = perturbRnd(Sigma, m) returns a m x nPars 
%        vector each of whose rows is a random sample from the perturbation 
%        distribution, whose parameters may depend on the covariance matrix
%        Sigma
%    evalDist -
%        function that when called as [d, output] = evalDist(Theta, threshold)
%        where Theta is a 1 x nPars vectors simulates the forward model under
%        parameters  Theta and returns the distance metric relative to the 
%        data (d) and a vector of model output (output). Including a threshold 
%        in the call allows optional early termination of simulations in
%        which the threshold is exceed 
%    hyperPar - structure of hyperparameters with the following fields:
%        nParticles - required number of samples from the posterior
%        acceptRate - acceptance probability for each SMC population
%        eps0 - initial tolerance for distance (used to implement early 
%               rejection - set to inf to accept fixed proportion in the 
%               rejection step regardless of distance)
%        epsTarg - target tolerance
%        fAcc - minimum MCMC acceptance rate
%        c - tuning parameter (there is a 1-c chance that all particles are
%            moved at least once); 0.01 in Carr et al
%        S0 - initial number of pilot MCMC steps
%    loadRejectStepFlag - if false algorithm will run from scratch; if true
%        the function will attempt to load the results of the rejection 
%        step from previously saved file in results/rejection_step.mat
%
% OUTPUTS
%    particles - table of particles sampled from the posterior with the
%    following fields:
%       Theta - nParticles x nPars matrix whose rows represent the accepted
%               particles
%       modelOutput - nParticles x k matrix whose rows represent the 
%                    corresponding model output that was returned for each particle
%       dist - nParticles x 1 vector of distances that were returned for 
%              each of the accepted particles


% This multiple of the identity matrix will be added to empirically
% calculated covariance matrices until all their eigenvalues are positive
small = 1e-8;


nRetain = ceil(hyperPar.acceptRate * hyperPar.nParticles);       % number of particles to retain at each stage


% Call model once with low early rejection tolerance to establish size of modelOutput
tmpTheta = priorRnd(1);
nPars = length(tmpTheta);
[~, tmpx] = evalDist(tmpTheta, 0.00001 * hyperPar.eps0);
[~, nOutputs] = size(tmpx);


if ~(loadRejectStepFlag)
    % Initial ABC rejection stage:
    fprintf('Running rejection step... \n');

    % Set up table of data for particles
    particles.Theta = zeros(hyperPar.nParticles, nPars);
    particles.dist = 2 * hyperPar.eps0 * ones(hyperPar.nParticles, 1);
    particles.modelOutput = zeros(hyperPar.nParticles, nOutputs);
    particles.ThetaPDF = nan(hyperPar.nParticles, 1);
    particles.nAccepted = zeros(hyperPar.nParticles, 1);
    particles.nSimulated = zeros(hyperPar.nParticles, 1);
    particles = struct2table(particles);


    % Iterate until at least required fraction of particles have met starting
    % threshold eps0
    acceptedFlag = zeros(hyperPar.nParticles, 1);
    iAttempts = 0;
    while sum(acceptedFlag) < nRetain

        % Pick out particles that have not yet been accepted and re-draw from
        % prior
        indToDo = find(~acceptedFlag);
        nToDo = length(indToDo);

        % Create temporary arrays with info for these re-drawn particles
        ThetaProp = priorRnd(nToDo);
        d = zeros(nToDo, 1);
        x = zeros(nToDo, nOutputs);
        iAttempts = iAttempts + nToDo;

        % Simulate model for these particles with early rejection tolerance eps0
        parfor iRep = 1:nToDo
          [d(iRep), x(iRep, :)] = evalDist(ThetaProp(iRep, :), hyperPar.eps0);
        end

        % Record results for the simulated particles in dist and modelOutput
        particles.Theta(indToDo, :) = ThetaProp;
        particles.dist(indToDo) = d;
        particles.modelOutput(indToDo, :) = x;

        % Check which particles have now met thresold
        acceptedFlag = (particles.dist < hyperPar.eps0);
        fprintf('    %i/%i (%.1f%%) met specified initial tolerance eps0 = %.3f\n', sum(acceptedFlag), iAttempts, 100*sum(acceptedFlag)/iAttempts, hyperPar.eps0)
    end

    % Sort particles by distance 
    [~, si] = sort(particles.dist);
    particles = particles(si, :);

    save('results/rejection_step.mat', 'particles');        % save results of rejection step to allow subsequent reloading if required
else
    fprintf('Loading rejection step... \n');
    load('results/rejection_step.mat', 'particles');        % load previously saved results of the rejection step
end

fprintf('...done.               dMax = %.3f\n', particles.dist(end));





% Series of SMC populations
pAcc = 1;               % make sure while loop runs at least once
St = hyperPar.S0;       % set initial number of MCMC trials
iPop = 0;               % population counter
while pAcc > hyperPar.fAcc & particles.dist(end) > hyperPar.epsTarg     % algorithm terminates when MCMC acceptance drops below fAcc OR all particles meet target distance epsTarg
   iPop = iPop+1;
   
    % Set new threshold according to number of accepted particles
    epst = particles.dist(nRetain);

    % Covariance matrix of accepted particles:
    Sigma = cov(particles.Theta(1:nRetain, :));
    while eigs(Sigma, 1, 'smallestreal') < small/2      % add small perturbation to ensure Sigma is positive definite
      Sigma = Sigma+small*eye(nPars); 
    end
    
   
   fprintf('SMC pop %i, epst = %.3f ', iPop, epst)
   % Replenish rejected particles by resampling from accepted particles
   nn = randi(nRetain, hyperPar.nParticles-nRetain, 1);
   particles(nRetain+1:end, :) = particles(nn, :);
   
   % Pre-evaluate prior PDF at each particle and store in table
   particles.ThetaLogPDF = priorLogPDF(particles.Theta);
   
   % Setup vectors to record the number of MCMC proposal points that are
   % simulated and the number that are accepted 
   particles.nAccepted = zeros(hyperPar.nParticles, 1);  
   particles.nSimulated = zeros(hyperPar.nParticles, 1);  

   % Propagate replenished particles through MCMC steps
   fprintf('MCMC steps (%i+', St)
   parfor iParticle = nRetain+1:hyperPar.nParticles
       particles(iParticle, :) = doMCMCsteps(particles(iParticle, :), priorLogPDF, perturbRnd, evalDist, Sigma, epst, St);
   end
   
   % Compute number of additional MCMC steps required
   pAcc = sum(particles.nAccepted)/(St*(hyperPar.nParticles-nRetain));
   Rt = ceil( log(hyperPar.c)/log(1-pAcc) );
   nExtra = max(0, Rt-St);
   
   % Run nExtra additional MCMC steps
   fprintf('%i=%i)... ', nExtra, St+nExtra)
   if nExtra > 0
       parfor iParticle = nRetain+1:hyperPar.nParticles
           particles(iParticle, :) = doMCMCsteps(particles(iParticle, :), priorLogPDF, perturbRnd, evalDist, Sigma, epst, nExtra);
       end
   end
   
   % Sort particles by distance 
    [~, si] = sort(particles.dist);
    particles = particles(si, :);
    
    % Calculate proportion of particles that are unique
    pUnique = size(unique(particles.Theta, 'rows'), 1)/hyperPar.nParticles;

    % Set number of pilot steps for next SMC population
    totAttempts = (St+nExtra)*(hyperPar.nParticles-nRetain);
    totSimulated = sum(particles.nSimulated);
    totAccepted = sum(particles.nAccepted);
    pAcc = totAccepted/totAttempts;
    Rt = ceil( log(hyperPar.c)/log(1-pAcc) );
    St = ceil(Rt/2);
    
    fprintf('done.  dMax = %.3f, %i/%i (%.3f%%) simulated, %i/%i (%.3f%%) accepted, pUnique = %.3f\n', particles.dist(end), totSimulated, totAttempts, 100*totSimulated/totAttempts, totAccepted, totAttempts, 100*pAcc, pUnique);
    
end




