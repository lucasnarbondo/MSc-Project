clear all
clc
WindCurtailed=[5.949 6.547 40.032 68.640 59.795 186.083 152.175 169.375 154.046 155.942 63.131 36.537];
WindProduction=[298.394 197.348 322.276 345.535 316.891 259.691 411.346 280.503 270.449 315.858 358.966 377.088];

WindPlot=[298.394, 5.949;
          197.348, 6.547;
          322.276, 40.032;
          345.535, 68.640;
          316.891, 59.795;
          259.691, 186.083;
          411.346, 152.175;
          280.503, 169.375;
          270.449, 154.046;
          315.858, 155.942;
          358.966, 63.131;
          377.088, 36.537];

figure
h=area(WindPlot);
h(1).FaceColor = [0 0.5 0.5];
h(2).FaceColor = [0 0.75 0.75];
grid
xticks([1 2 3 4 5 6 7 8 9 10 11 12])
xlim([1 12])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel('Energy (GWh)','Color','k')
legend('Wind Energy Produced','Wind Energy Curtailed','Location','northwest','Interpreter','latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network\Images';
baseFileName =  sprintf("WindCurtailed2017.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);