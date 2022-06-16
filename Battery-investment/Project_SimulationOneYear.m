%clear all
%clc

%load GurobiDataYear

%%%%%%%%%%%%%%%%%%%%
%% BATTERIES DATA %%
%%%%%%%%%%%%%%%%%%%%

Bat_investment_cost=90; %USD/kWh
Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000; %USD/MWh/year

Maintenance_costs=0.05*Bat_ANNUITIZED; %Maintenance assumed to be 5% of the capital costs

Bat_rate=150;   %Charge/Discharge rate
Bat_ef=0.94;   %Eficiency

%%%%%%%%%%%%%%%%%%%%%
%% LOAD SPOT PRICE %%
%%%%%%%%%%%%%%%%%%%%%

Year=2016;
FileName= ['Spot' num2str(Year)];
folder='C:\Users\lucas\Documents\Imperial College\Project\Battery investment with EFR\Spot';
fullFileName = fullfile(folder, FileName);

SpotPrice=xlsread(fullFileName,'B1:Y365');
SpotVector=reshape(transpose(SpotPrice),[],1); %Re-arrange data as vector

%%%%%%%%%%%%%%%%
%% OTHER DATA %%
%%%%%%%%%%%%%%%%

%Transmission_fees=105.38/35*1000*12; %original fees from MIEM (/12 per month test)
Transmission_fees=0; %original fees from MIEM (/12 per month test)

TotalSimulation=length(SpotVector);   %simlulate over a year
aux_zeros=zeros(1,TotalSimulation);   %Auxiliar Vector for algorithm

%%%%%%%%%%%%%
%% PROBLEM %%
%%%%%%%%%%%%%

%% Define Objective Function

model.obj=[(-Bat_ANNUITIZED-Maintenance_costs) aux_zeros -SpotVector' SpotVector']; %Maximize revenue (Profits-Costs-Investment)
model.objcon=-Transmission_fees*Bat_rate; %Transmission fees as a constant

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
BatterySize=result.x(1);
StorageLevel=result.x(2:1+TotalSimulation);
ChargeRate=result.x(1+TotalSimulation:1+2*TotalSimulation);
DisChargeRate=result.x(1+2*TotalSimulation:1+3*TotalSimulation);

%% Plots
%Plots
% figure
% plot(StorageLevel,'linewidth',0.3)
% ylabel('StorageLevel (MWh)','Color','k')
% xticks([0 720  2160 2880 3600 4320 5040 5760 6480 7200 37920])
% xticklabels({'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
% pbaspect([4 1 1])
% ax = gca;
% ax.FontSize = 12;
% ax.FontName = 'Times';
% set(gcf, 'Color', 'w');

fprintf('Optimal Battery Capacity = %.2f \n',result.x(1))