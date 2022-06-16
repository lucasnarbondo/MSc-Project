% LOAD WIND AND DEMAND DATA
% Information from a whole week
% Time Period 10 minutes

%Case 2

clear all
clc

SG=xlsread('GenerationDemandData_Dry_2017','B3460:B4467');
Bonete=xlsread('GenerationDemandData_Dry_2017','C3460:C4467');
Baygorria=xlsread('GenerationDemandData_Dry_2017','D3460:D4467');
Palmar=xlsread('GenerationDemandData_Dry_2017','E3460:E4467');
RioNegro=Bonete+Baygorria+Palmar;
Wind=xlsread('GenerationDemandData_Dry_2017','F3460:F4467');
Demand=xlsread('GenerationDemandData_Dry_2017','M3460:M4467');
Solar=xlsread('GenerationDemandData_Dry_2017','G3460:G4467');
Biomass=xlsread('GenerationDemandData_Dry_2017','I3460:I4467');
ExpArg=xlsread('GenerationDemandData_Dry_2017','N3460:N4467');
ExpBrasil=xlsread('GenerationDemandData_Dry_2017','O3460:O4467');
Thermal=xlsread('GenerationDemandData_Dry_2017','H3460:H4467');

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
% baseFileName =  sprintf("ThermalGenerationDry2017-Case2.pdf");
% fullFileName = fullfile(folder, baseFileName);
% export_fig(fullFileName);

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
baseFileName =  sprintf("ThermalHistogramDry2017-Case2.pdf");
fullFileName = fullfile(folder, baseFileName);
%export_fig(fullFileName);

save('GenerationAndDemandData_Dry_2017','SaltoGrandeAvg','RioNegroAvg','WindAvg','DemandAvg','SolarAvg','BiomassAvg','ExpArgAvg','ExpBrasilAvg','ThermalAvg')
