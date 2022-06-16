clear all
clc
ThermalEnergy=[21.4 26.3 4.1 1.2 5.5 45.7 7.4 6.7 35.1 2.2+9.2 11.8 44.1];
bar(ThermalEnergy)
grid
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Energy (GWh)','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network\Images';
baseFileName =  sprintf("ThermalGeneration2017.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);