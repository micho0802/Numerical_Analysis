clear all, clc
A = readtable('US-COVID-11-7-22.xlsx');
A1 = A(1:250,"activeCasesSinceFeb15_2020");
A2 = A(251:500,"activeCasesSinceFeb15_2020");
A3 = A(501:750,"activeCasesSinceFeb15_2020");
A4 = A(751:996,"activeCasesSinceFeb15_2020");
X = 1:1200;

%Custom function

y = @(x) 1.497e+07*exp(-((x-715.8 )/25.28).^2) + 3.565e+06 *exp(-((x-793.8)/204.1).^2) - 2.479e+06*exp(-((x-793.6)/34.39).^2) + 6.48e+06*exp(-((x-335.7)/57.13).^2) + 1.876e+06*exp(-((x-201.1)/101.6).^2) + 3.846e+06*exp(-((x-579.9)/40.19).^2) + 1.864e+06*exp(-((x-441)/39.51).^2) + 1.419e+06*exp(-((x-909.8)/38.69).^2);
fplot(y,[1,1200]);