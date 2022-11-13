function y = normLogPDF(x, Mu, Sigma)

% normLogPDF - log of normal PDF
%   y = normLogPDF(x, Mu, Sigma) returns the log of the PDF of the univariate
%   normal distribution with mean Mu and variance Sigma evaluated at each
%   value in the array x

y = -0.5*((x - Mu)./Sigma).^2 - 0.5*log(2) - 0.5*log(pi) - log(Sigma); 


