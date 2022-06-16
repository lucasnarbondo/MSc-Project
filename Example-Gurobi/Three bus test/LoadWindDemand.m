% LOAD WIND AND DEMAND DATA
% Information from a whole month
% Time Period 10 minutes
% 
% File WindDemandData.xls

clear all
clc

WindData=xlsread('WindDemandData','F4:F4467');
DemandData=xlsread('WindDemandData','M4:M4467');

TimeStep=10/60;  %Data every 10 minutes
TotalSimulation=24/TimeStep;
NumberDays=31;

WindAvg=zeros(1,TotalSimulation);
DemandAvg=zeros(1,TotalSimulation);


for k=1:length(WindAvg)
    aux_wind=WindData(k);
    aux_demand=DemandData(k);
    for i=1:NumberDays-1
        aux_wind=aux_wind+WindData(i*TotalSimulation+k);
        aux_demand=aux_demand+DemandData(i*TotalSimulation+k);
    end
    WindAvg(k)=aux_wind/NumberDays;
    DemandAvg(k)=aux_demand/NumberDays;
end

save('WindAndDemandData','WindAvg','DemandAvg')
