clear all
clc
ThermalEnergy=[71.1 12.6 58.5 13.6 4.2 59.4 68.04 23.9 16.7 14.1 2.3 5.3];
bar(ThermalEnergy)
grid
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Energy (GWh)','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network\Images';
baseFileName =  sprintf("Thermal Generation.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);