figure
%plot(Price,BatterySize2018,'k','linewidth',0.3);
plot(Price,BatterySize2018,'k','linewidth',0.3);
hold on
plot(Price,BatterySize2017,'color',[0 0.4470 0.7410],'linewidth',0.3);
plot(Price,BatterySize2016,'color',[0.9290 0.6940 0.1250],'linewidth',0.3);
plot(Price,BatterySize2015,'color',[0.4940 0.1840 0.5560],'linewidth',0.3);
plot(Price,BatterySize2014,'color',[0.4660 0.6740 0.1880],'linewidth',0.3);
plot(Price,BatterySize2013,'color',[0.8500 0.3250 0.0980],'linewidth',0.3);
hold off
xlim([80 230])
ylim([0 800]);
grid
xlabel('Price (USD/kWh)','Color','k') 
ylabel('Capacity (MWh)','Color','k')
title('Battery Size as a function of price comparison','Color','k')
legend('2018','2017','2016','2015','2014','2013','Interpreter','latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Battery investment\Images';
baseFileName =  sprintf("Comparison Battery Size All.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);