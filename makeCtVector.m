function Ct = makeCtVector(Theta, par) 

nTheta = length(Theta);

Ct = ones(1, par.tEnd+1);
 

Ct(1:8) = Theta(1);
for ii = 2:nTheta
    Ct(7*ii-5:7*ii+1) = Theta(ii); 
end
Ct(7*nTheta+2:end) = Theta(end);

