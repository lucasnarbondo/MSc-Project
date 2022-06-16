% LOAD WIND AND DEMAND DATA
% Information from a whole week
% Time Period 10 minutes

%Case 2

clear all
clc

SG=xlsread('GenerationDemandData_Wet_2017-2','B436:B1587');
Bonete=xlsread('GenerationDemandData_Wet_2017-2','C436:C1587');
Baygorria=xlsread('GenerationDemandData_Wet_2017-2','D436:D1587');
Palmar=xlsread('GenerationDemandData_Wet_2017-2','E436:E1587');
RioNegro=Bonete+Baygorria+Palmar;
Wind=xlsread('GenerationDemandData_Wet_2017-2','F436:F1587');
Demand=xlsread('GenerationDemandData_Wet_2017-2','M436:M1587');
Solar=xlsread('GenerationDemandData_Wet_2017-2','G436:G1587');
Biomass=xlsread('GenerationDemandData_Wet_2017-2','I436:I1587');
ExpArg=xlsread('GenerationDemandData_Wet_2017-2','N436:N1587');
ExpBrasil=xlsread('GenerationDemandData_Wet_2017-2','O436:O1587');
Thermal=xlsread('GenerationDemandData_Wet_2017-2','H436:H1587');

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

TotalDemand=DemandAvg+ExpArgAvg+ExpBrasilAvg;
TotalPower=WindAvg+SolarAvg+BiomassAvg+SaltoGrandeAvg+RioNegroAvg;

for k=1:length(TotalDemand)
    if TotalDemand(k)<TotalPower(k)
        DemandAvg(k)=WindAvg(k)+SolarAvg(k)+BiomassAvg(k)+SaltoGrandeAvg(k)+RioNegroAvg(k)-ExpArgAvg(k)-ExpBrasilAvg(k);
    end
end

%Plot Thermal
figure
plot(ThermalAvg,'k','linewidth',0.3)
hold on
xticks([12 36 60 84 108 132 156])
xticklabels({'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7'})
xtickangle(45)
ylabel('Energy (MWh)','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network Investment\Images';
baseFileName =  sprintf("ThermalGenerationWet2017-Case2.pdf");
fullFileName = fullfile(folder, baseFileName);
%export_fig(fullFileName);

figure
histogram(Thermal,25,'facecolor','k');
grid
xlabel('MW','Color','k') 
ylabel('Number of Times','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network Investment\Images';
baseFileName =  sprintf("ThermalHistogramWet2017-Case2.pdf");
fullFileName = fullfile(folder, baseFileName);
% export_fig(fullFileName);

save('GenerationAndDemandData_Wet_2017_2','SaltoGrandeAvg','RioNegroAvg','WindAvg','DemandAvg','SolarAvg','BiomassAvg','ExpArgAvg','ExpBrasilAvg','ThermalAvg')
