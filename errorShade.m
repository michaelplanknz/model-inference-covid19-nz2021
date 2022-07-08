function errorShade(x, y, ci, clr, ls, mkr)

x = reshape(x, 1, length(x));
y = reshape(y, 1, length(y));
[n, m] = size(ci);
if n == length(x) & m ~= length(x)
    ci = ci.';
end
[n, m] = size(ci);
nBands = n/2;
if nBands > 7
    fprintf(    'WARNING: errorShade cannot plot >7 bands\n\n')
end

clrCoeff = 0.75 - 0.2*(nBands/7)*linspace(-1, 1, nBands);

for iBand = 1:nBands
    ind1 = iBand;
    ind2 = 2*nBands+1-iBand;
    inLowerFlag = ~isnan(ci(ind1, :));
    inUpperFlag = ~isnan(ci(ind2, :));

    xShade = [x(inLowerFlag), fliplr(x(inUpperFlag))];
    yShade = [ci(ind1, inLowerFlag), fliplr(ci(ind2, inUpperFlag))];
    fill(xShade, yShade, clrCoeff(iBand) + (1-clrCoeff(iBand))*clr, 'LineStyle', 'none', 'HandleVisibility', 'off', 'FaceAlpha', 0.5)
    hold on
end
plot(x, y, 'LineStyle', ls, 'Marker', mkr, 'Color', clr)
