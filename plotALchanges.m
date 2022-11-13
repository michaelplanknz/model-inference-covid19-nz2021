function plotALchanges()

xline(datetime('18AUG2021', 'InputFormat', 'ddMMMyyyy'), 'k--', 'AL4');
xline(datetime('22SEP2021', 'InputFormat', 'ddMMMyyyy'), 'k--', 'AL3');
xline(datetime('06OCT2021', 'InputFormat', 'ddMMMyyyy'), 'k--', 'AL3/S1');
xline(datetime('10NOV2021', 'InputFormat', 'ddMMMyyyy'), 'k--', 'AL3/S2');
xline(datetime('03DEC2021', 'InputFormat', 'ddMMMyyyy'), 'k--', 'CPF red');
xline(datetime('31DEC2021', 'InputFormat', 'ddMMMyyyy'), 'k--', 'CPF orange');
