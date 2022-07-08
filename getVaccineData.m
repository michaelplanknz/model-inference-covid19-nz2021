function [dates, V1, V2] = getVaccineData(vaxDataDate)

fLbl = datestr(vaxDataDate, 'yyyy-mm-dd');

fNameDates = "data/dates_vaccine.csv";
fName1 = "data/firstDoseProp_akl_" + fLbl + ".csv";
fName2 = "data/secondDoseProp_akl_" + fLbl + ".csv";

fprintf('   Loading vaccination data:    %s\n                                %s\n                                %s\n', fNameDates, fName1, fName2)
dates = readmatrix(fNameDates, 'OutputType', 'datetime');
V1 = readmatrix(fName1).';
V2 = readmatrix(fName2).';



