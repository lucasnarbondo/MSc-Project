% clear all
% clc
% 
% load GurobiDataYear

%%%%%%%%%%%%%%%%%%%%%
%% LOAD SPOT PRICE %%
%%%%%%%%%%%%%%%%%%%%%

Year=2013;
FileName= ['Spot' num2str(Year)];
folder='C:\Users\lucas\Documents\Imperial College\Project\Battery investment\Spot';
fullFileName = fullfile(folder, FileName);
SpotPrice=xlsread(fullFileName,'B1:Y365');
SpotVector=reshape(transpose(SpotPrice),[],1); %Re-arrange data as vector

%%%%%%%%%%%%%%%%
%% OTHER DATA %%
%%%%%%%%%%%%%%%%

Transmission_fees=105.38/35*1000*12;
TotalSimulation=length(SpotVector);   %simlulate over a year
aux_zeros=zeros(1,TotalSimulation);   %Auxiliar Vector for algorithm

%%%%%%%%%%%%%%%%%%%%
%% BATTERIES DATA %%
%%%%%%%%%%%%%%%%%%%%

Bat_investment_cost=80; %USD/kWh
Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value
Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000; %USD/MWh/year

%%%%%%%%%%%%%
%% PROBLEM %%
%%%%%%%%%%%%%

BatterySize=zeros(1,100);
Rate=zeros(1,100);

for k=1:100

Bat_rate=50+2*(k-1); %USD/MWh/year;   %Charge/Discharge rate 
Rate(k)=Bat_rate;

%% Define Objective Function

model.obj=[(-Bat_ANNUITIZED-Maintenance_costs) aux_zeros -SpotVector' SpotVector']; %Maximize revenue (Profits-Costs-Investment)
model.objcon=-Transmission_fees*Bat_rate;

%% Bounded Equations

LowerBound=[0 aux_zeros aux_zeros aux_zeros];

aux_inf=inf*ones(1,TotalSimulation);
aux_rate=Bat_rate*ones(1,TotalSimulation);

UpperBound=[inf aux_inf aux_rate aux_rate];
    
model.lb=LowerBound';
model.ub=UpperBound';

%% Gurobi Optimization
tic
result=gurobi(model,params);
toc

%% Results
BatterySize2013(k)=result.x(1);
end


%% Plots
% figure
% plot(Price,BatterySize,'k','linewidth',0.3);
% xlim([80 230])
% ylim([0 800]);
% grid
% xlabel('Price (USD/kWh)','Color','k') 
% ylabel('Capacity (MWh)','Color','k')
% title(['Battery Size as a function of price -  ' num2str(Year)],'Color','k')
% ax = gca;
% ax.FontSize = 12;
% ax.FontName = 'Times';
% set(gcf, 'Color', 'w');
% folder = 'C:\Users\lucas\Documents\Imperial College\Project\Battery investment\Images';
% baseFileName =  sprintf("Bat Size Spot %d.pdf",Year);
% fullFileName = fullfile(folder, baseFileName);
% export_fig(fullFileName);
