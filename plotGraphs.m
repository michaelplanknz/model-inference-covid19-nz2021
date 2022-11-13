function plotGraphs(epiData, fittedToDate)


fName = sprintf('results/results_fit_to_%s.mat', datestr(fittedToDate));
load(fName); 

Reff = Rvt.*Ct.*(1-TTIQeff_time);

epiNow = readtable('data/EpiNow2_estimates_AllRegions_MOHdata_20220126.csv');
epiNowFull = epiNow(epiNow.variable == "R" & epiNow.type == "estimate" , :);



qt = [0.1 0.2 0.3 0.4 0.6 0.7 0.8 0.9];      % quantiles to plot as uncertainty bands
clr = [0 0 1];
ls = '-';
mkr = 'none';
indPlot = 5:7:length(t);

t = datetime(t, 'ConvertFrom', 'datenum');

h = figure(1);
h.Position = [    291         111        1141         885];
subplot(3, 2, 1)
errorShade(t(indPlot), median(Ct(:, indPlot)), quantile(Ct(:, indPlot), qt), clr, ls, mkr);
plotALchanges();
ylabel('control effect C(t)')
xlim([t(1),t(end)])
ylim([0 1.2])
title('(a)')
subplot(3, 2, 2)
errorShade(t(indPlot), median(Reff(:, indPlot)), quantile(Reff(:, indPlot), qt), clr, ls, mkr);
hold on
plot(epiNowFull.date.', epiNowFull.median.', 'r-', epiNowFull.date.', epiNowFull.upper_90.', 'r--', epiNowFull.date.', epiNowFull.lower_90.', 'r--');
plotALchanges();
yline(1, 'k:');
ylabel('Reff')
ylim([0 2])
xlim([t(1),t(end)])
ylim([0 2])
title('(b)')

subplot(3, 2, 3)
errorShade(t, median(nIsolAll), quantile(nIsolAll, qt), clr, ls, mkr)
hold on
fittedFlag = datenum(epiData.Date) <= fittedToDate;
plot(epiData.Date(fittedFlag), epiData.nCases(fittedFlag), 'go', epiData.Date(~fittedFlag), epiData.nCases(~fittedFlag), 'bo')
plotALchanges();
ylabel('daily reported cases')
xlim([t(1),t(end)])
title('(c)')
subplot(3, 2, 4)
errorShade(t, median(nHospAll), max(1, quantile(nHospAll, qt)), clr, ls, mkr)
hold on
plotALchanges();
ylabel('new daily hospital admissions')
xlim([t(1),t(end)])
title('(d)')
subplot(3, 2, 5)
errorShade(t, median(cumsum(nDeathsAll, 2)), max(1, quantile(cumsum(nDeathsAll, 2), qt)), clr, ls, mkr)
hold on
plotALchanges();
ylabel('cumulative deaths')
xlim([t(1),t(end)])
title('(e)')
subplot(3, 2, 6)
errorShade(t, median(nBedsAll), quantile(nBedsAll, qt ), clr, ls, mkr)
hold on
plot(epiData.Date, epiData.hospBeds, 'bo')
plotALchanges();
ylabel('hospital beds occupied')
xlim([t(1),t(end)])
title('(f)')

saveas(h, sprintf('figures/results_fit_to_%s.png', datestr(fittedToDate)));









































