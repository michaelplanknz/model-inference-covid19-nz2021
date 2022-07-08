function vaccDoses = genSeedVaxStatus_copyPopn(cases, nSeedCases, par)

pv1 = ((1-par.VEi1)*par.cov1(1, cases.ageGroup(1:nSeedCases)) ./ (1-par.VEi1*par.cov1(1, cases.ageGroup(1:nSeedCases))-par.VEi2*par.cov2(1, cases.ageGroup(1:nSeedCases)))).';
pv2 = ((1-par.VEi2)*par.cov2(1, cases.ageGroup(1:nSeedCases)) ./ (1-par.VEi1*par.cov1(1, cases.ageGroup(1:nSeedCases))-par.VEi2*par.cov2(1, cases.ageGroup(1:nSeedCases)))).';
r = rand(nSeedCases, 1);
vaccDoses = (r < pv2) + (r < pv1+pv2);
