
function pertSigma = getPertSigma(priorSigma, popnSigma)

% Calculate covariance matrix to define the pertrubation distribution   

Alpha = 0; Beta = 2;

pertSigma = zeros(size(priorSigma));
pertSigma(1, 1) = Alpha*priorSigma(1, 1) + Beta*popnSigma(1, 1);
pertSigma(2:end, 2:end) = (Alpha + Beta*trace(popnSigma(2:end, 2:end))/trace(priorSigma(2:end, 2:end)) ) * priorSigma(2:end, 2:end);


