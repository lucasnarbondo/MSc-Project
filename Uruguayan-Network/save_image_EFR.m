figure
plot(Price,BaterySize_EFR,'k','linewidth',0.3);
hold on
plot(Price,BatterySize,'k:','linewidth',0.3);
hold off
xlim([80 180])
ylim([0 600]);
grid
xlabel('Price (USD/kWh)','Color','k') 
ylabel('Capacity (MWh)','Color','k')
title('Battery Size as a function of price comparison','Color','k')
legend('With EFR','Without EFR','Interpreter','latex')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Battery investment with EFR\Images';
baseFileName =  sprintf("Comparison Battery EFR.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);