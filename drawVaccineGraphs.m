clear all
close all


vaxDataDate = [datenum('14OCT2021'); datenum('25NOV2021'); datenum('17JAN2022') ];

par1 = getPar(vaxDataDate(1));
par2 = getPar(vaxDataDate(2));
par3 = getPar(vaxDataDate(3));

[nVaxDates, ~] = size(par3.cov1);

t = par1.date0:par1.date0+nVaxDates-1;
elig = [0 0 3/5 ones(1, 13)];
cov1_12plus_1 = sum(par1.cov1.*par1.popCount.', 2)/( sum(elig.*par1.popCount.') );
cov2_12plus_1 = sum(par1.cov2.*par1.popCount.', 2)/( sum(elig.*par1.popCount.') );
cov1_12plus_2 = sum(par2.cov1.*par1.popCount.', 2)/( sum(elig.*par1.popCount.') );
cov2_12plus_2 = sum(par2.cov2.*par1.popCount.', 2)/( sum(elig.*par1.popCount.') );
cov1_12plus_3 = sum(par3.cov1.*par1.popCount.', 2)/( sum(elig.*par1.popCount.') );
cov2_12plus_3 = sum(par3.cov2.*par1.popCount.', 2)/( sum(elig.*par1.popCount.') );



% get dominant eigenvalue of vaccinated NGM through time to calculater
% effect of vaccination on reproduciton number
Rvi = zeros(1, length(t));
for iDate = 1:length(t)
    Rvi_1(iDate) = calcVaxedR(1-par1.cov1(iDate, :)-par1.cov2(iDate, :), par1.cov1(iDate, :), par1.cov2(iDate, :), par1);
    Rvi_2(iDate) = calcVaxedR(1-par2.cov1(iDate, :)-par2.cov2(iDate, :), par2.cov1(iDate, :), par2.cov2(iDate, :), par2);
    Rvi_3(iDate) = calcVaxedR(1-par3.cov1(iDate, :)-par3.cov2(iDate, :), par3.cov1(iDate, :), par3.cov2(iDate, :), par3);
end


figure
subplot(1, 2, 1)
plot(t, cov1_12plus_1, 'b-.', t, cov2_12plus_1, 'r-.', t, cov1_12plus_1+cov2_12plus_1, 'k-.')
hold on
plot(t, cov1_12plus_2, 'b--', t, cov2_12plus_2, 'r--', t, cov1_12plus_2+cov2_12plus_2, 'k--')
plot(t, cov1_12plus_3, 'b-', t, cov2_12plus_3, 'r-', t, cov1_12plus_3+cov2_12plus_3, 'k-')
ylabel('proportion of 12+')
xline(vaxDataDate(1), 'k:')
xline(vaxDataDate(2), 'k:')
xline(vaxDataDate(3), 'k:')
ylim([0 1])
xlim([t(1) t(1)+par1.tEnd])
datetick('x', 'dd-mmm', 'keeplimits')
title('(a)')
h = gca;
legend([h.Children(6), h.Children(5), h.Children(4)], "1 dose", "2 doses", "at least 1 dose");
subplot(1, 2, 2)
plot(t, Rvi_1/par1.R0, 'b-.')
hold on
plot(t, Rvi_2/par1.R0, 'b--')
plot(t, Rvi_3/par1.R0, 'b-')
xline(vaxDataDate(1), 'k:')
xline(vaxDataDate(2), 'k:')
xline(vaxDataDate(3), 'k:')
ylim([ 0 1])
xlim([t(1) t(1)+par1.tEnd])
datetick('x', 'dd-mmm', 'keeplimits')
ylabel('relative reproduction number')
title('(b)')
legend('Oct 2021 data', 'Nov 2021 data', 'Jan 2022 data')

