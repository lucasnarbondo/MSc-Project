% LOAD WIND AND DEMAND DATA
% Information from a whole month
% Time Period 10 minutes
% 
% File 2013.xls

clear all
clc

WindData=xlsread('GenerationDemandData','F4:F4467');
DemandData=xlsread('GenerationDemandData','M4:M4467');
SolarData=xlsread('GenerationDemandData','G4:G4467');
BiomassData=xlsread('GenerationDemandData','I4:I4467');
ExpArgData=xlsread('GenerationDemandData','N4:N4467');
ExpBrasilData=xlsread('GenerationDemandData','O4:O4467');
ThermalData=xlsread('GenerationDemandData','H4:H4467');

TimeStep=10/60;  %Data every 10 minutes
TotalSimulation=24/TimeStep;
NumberDays=31;

WindAvg=zeros(1,TotalSimulation);
DemandAvg=zeros(1,TotalSimulation);
SolarAvg=zeros(1,TotalSimulation);
BiomassAvg=zeros(1,TotalSimulation);
ExpArgAvg=zeros(1,TotalSimulation);
ExpBrasilAvg=zeros(1,TotalSimulation);
ThermalAvg=zeros(1,TotalSimulation);

for k=1:length(WindAvg)
    aux_wind=WindData(k);
    aux_demand=DemandData(k);
    aux_solar=SolarData(k);
    aux_biomass=BiomassData(k);
    aux_exparg=ExpArgData(k);
    aux_expbrasil=ExpBrasilData(k);
    aux_thermal=ExpBrasilData(k);
    for i=1:NumberDays-1
        aux_wind=aux_wind+WindData(i*TotalSimulation+k);
        aux_demand=aux_demand+DemandData(i*TotalSimulation+k);
        aux_solar=aux_solar+SolarData(i*TotalSimulation+k);
        aux_biomass=aux_biomass+BiomassData(i*TotalSimulation+k);
        aux_exparg=aux_exparg+ExpArgData(i*TotalSimulation+k);
        aux_expbrasil=aux_expbrasil+ExpBrasilData(i*TotalSimulation+k);
        aux_thermal=aux_thermal+ThermalData(i*TotalSimulation+k);
    end
    WindAvg(k)=aux_wind/NumberDays;
    DemandAvg(k)=aux_demand/NumberDays;
    SolarAvg(k)=aux_solar/NumberDays;
    BiomassAvg(k)=aux_biomass/NumberDays;
    ExpArgAvg(k)=aux_exparg/NumberDays;
    ExpBrasilAvg(k)=aux_expbrasil/NumberDays;
    ThermalAvg(k)=aux_thermal/NumberDays;
end

save('GenerationAndDemandData','WindAvg','DemandAvg','SolarAvg','BiomassAvg','ExpArgAvg','ExpBrasilAvg','ThermalAvg')