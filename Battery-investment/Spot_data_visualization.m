%%%%%%%%%%%%%%%%%%%%%%%%
%% DATA VISUALIZATION %%
%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%%%%%%%%%%%%%%%%%%%%%
%% LOAD SPOT PRICE %%
%%%%%%%%%%%%%%%%%%%%%

Year=2013;
FileName= ['Spot' num2str(Year)];
folder='C:\Users\lucas\Documents\Imperial College\Project\Battery investment with EFR\Spot';
fullFileName = fullfile(folder, FileName);

SpotPrice=xlsread(fullFileName,'B1:Y365');
SpotPrice=transpose(SpotPrice);

%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZATION PRICE %%
%%%%%%%%%%%%%%%%%%%%%%%%%

c=gray;
c=flipud(c);
figure
imagesc(SpotPrice)
colormap(c)
c=colorbar;
c.Label.String='USD/MWh';
yticks([8 16])
yticklabels({'8:00','16:00'})
xticks([0 31 59 90 120 151 181 212 243 273 304 334])
xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
pbaspect([4 1 1])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Battery investment\Images';
baseFileName =  sprintf("SpotVisualization%d.pdf",Year);
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALIZATION DIFFERENCE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Difference=max(SpotPrice)-min(SpotPrice);
figure
histogram(Difference,25,'facecolor','k');
grid
xlabel('USD/MWh','Color','k') 
ylabel('Number of Days','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Battery investment\Images';
baseFileName =  sprintf("SpotDifference%d.pdf",Year);
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);