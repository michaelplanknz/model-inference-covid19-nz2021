function [V1out, V2out] = convertFirstDoses(dates, V1, V2, startDate, tChar, par)

unconverted = V1(end, :);
conversions = (dates' >= startDate) .* min(1, (dates'-startDate)/(tChar)) .* unconverted;
fprintf('   Converting first doses to second in the %i days from %s\n', tChar, datestr(startDate))
%fprintf('       age group %i: %.2f%%\n', [1:par.nAgeGroups; 100*unconverted])

V1out = V1-conversions;
V2out = V2+conversions;

