function Rvi = calcVaxedR(susFrac0, susFrac1, susFrac2, par)    

infBlock = susFrac0 + (1-par.VEi1)*susFrac1 + (1-par.VEi2)*susFrac2;    % row
transBlockNum0 =                               (par.IDR + par.cSub*(1-par.IDR))                           .* susFrac0';
transBlockNum1 = (1-par.VEi1) * (1-par.VEt1) * ((1-par.VEs1)*par.IDR + par.cSub*(1-(1-par.VEs1)*par.IDR)) .* susFrac1';
transBlockNum2 = (1-par.VEi2) * (1-par.VEt2) * ((1-par.VEs2)*par.IDR + par.cSub*(1-(1-par.VEs2)*par.IDR)) .* susFrac2';  
transBlock = (transBlockNum0+transBlockNum1+transBlockNum2)./infBlock';                                     % column
NGMv = par.NGMclin .* infBlock .* transBlock;       % NGM with vaccine coverage on specified date
Rvi = eigs(NGMv, 1);    

