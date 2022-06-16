% LOAD WIND AND DEMAND DATA
% Information from a whole week
% Time Period 10 minutes

clear all
clc

SG=xlsread('GenerationDemandData2','B4:B1011');
Bonete=xlsread('GenerationDemandData2','C4:C1011');
Baygorria=xlsread('GenerationDemandData2','D4:D1011');
Palmar=xlsread('GenerationDemandData2','E4:E1011');
RioNegro=Bonete+Baygorria+Palmar;
Wind=xlsread('GenerationDemandData2','F4:F1011');
Demand=xlsread('GenerationDemandData2','M4:M1011');
Solar=xlsread('GenerationDemandData2','G4:G1011');
Biomass=xlsread('GenerationDemandData2','I4:I1011');
ExpArg=xlsread('GenerationDemandData2','N4:N1011');
ExpBrasil=xlsread('GenerationDemandData2','O4:O1011');
Thermal=xlsread('GenerationDemandData2','H4:H1011');

%Do the 10-min average for each hour
S  = numel(Wind);

xx = reshape(SG(1:S - mod(S, 6)), 6, []);
SaltoGrandeAvg = sum(xx, 1).' / 6;

xx = reshape(RioNegro(1:S - mod(S, 6)), 6, []);
RioNegroAvg = sum(xx, 1).' / 6;

xx = reshape(Wind(1:S - mod(S, 6)), 6, []);
WindAvg = sum(xx, 1).' / 6;

xx = reshape(Demand(1:S - mod(S, 6)), 6, []);
DemandAvg = sum(xx, 1).' / 6;

xx = reshape(Solar(1:S - mod(S, 6)), 6, []);
SolarAvg = sum(xx, 1).' / 6;

xx = reshape(Biomass(1:S - mod(S, 6)), 6, []);
BiomassAvg = sum(xx, 1).' / 6;

xx = reshape(ExpArg(1:S - mod(S, 6)), 6, []);
ExpArgAvg = sum(xx, 1).' / 6;

xx = reshape(ExpBrasil(1:S - mod(S, 6)), 6, []);
ExpBrasilAvg = sum(xx, 1).' / 6;

xx = reshape(Thermal(1:S - mod(S, 6)), 6, []);
ThermalAvg = sum(xx, 1).' / 6;


save('GenerationAndDemandData2','SaltoGrandeAvg','RioNegroAvg','WindAvg','DemandAvg','SolarAvg','BiomassAvg','ExpArgAvg','ExpBrasilAvg','ThermalAvg')
