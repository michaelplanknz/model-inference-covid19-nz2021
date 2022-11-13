function particle = doMCMCsteps(particle, priorLogPDF, perturbRnd, evalDist, Sigma, epst, nSteps)

for kStep = 1:nSteps
    ThetaProp = particle.Theta + perturbRnd(Sigma, 1);
    ThetaPropLogPDF = priorLogPDF(ThetaProp);
    MHratio = exp(ThetaPropLogPDF - particle.ThetaLogPDF  );
    if rand < MHratio                          % Don't need to simulate model if MH aapetence is 'no'
      [d, x] = evalDist(ThetaProp, epst);
      particle.nSimulated = particle.nSimulated+1;
      if d < epst 
          particle.Theta = ThetaProp;
          particle.modelOutput = x;
          particle.dist = d;
          particle.ThetaLogPDF = ThetaPropLogPDF;
          particle.nAccepted = particle.nAccepted+1;
      end
    end
end
