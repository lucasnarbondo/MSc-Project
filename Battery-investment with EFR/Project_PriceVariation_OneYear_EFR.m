%clear all
%clc

%load GurobiDataYear_EFR

%%%%%%%%%%%%%%%%%%%%%
%% LOAD SPOT PRICE %%
%%%%%%%%%%%%%%%%%%%%%

Year=2017;
FileName= ['Spot' num2str(Year)];
folder='C:\Users\lucas\Documents\Imperial College\Project\Battery investment with EFR\Spot';
fullFileName = fullfile(folder, FileName);

SpotPrice=xlsread(fullFileName,'B1:Y365');
SpotVector=reshape(transpose(SpotPrice),[],1); %Re-arrange data as vector

%%%%%%%%%
%% EFR %%
%%%%%%%%%

Time_response=15/60; %15 minutes supplying response
Payment_EFR=12; %Payment for availability to provide EFR in USD/MW/h

%%%%%%%%%%%%%%%%
%% OTHER DATA %%
%%%%%%%%%%%%%%%%

%Transmission_fees=105.38/35*1000*12; %original fees from MIEM
Transmission_fees=0;

TotalSimulation=length(SpotVector);   %simlulate over a year
aux_zeros=zeros(1,TotalSimulation);   %Auxiliar Vector for algorithm

%%%%%%%%%%%%%%%%%%%%
%% BATTERIES DATA %%
%%%%%%%%%%%%%%%%%%%%

Bat_rate=100;   %Charge/Discharge rate

Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

%%%%%%%%%%%%%
%% PROBLEM %%
%%%%%%%%%%%%%

BatterySize2017_EFR=zeros(1,100);
Price=zeros(1,100);

for k=1:100

Bat_investment_cost=40+1.5*(k-1);
Price(k)=Bat_investment_cost;

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000; %USD/MWh/year

Maintenance_costs=0.05*Bat_ANNUITIZED;

%% Define Objective Function

model.obj=[(-Bat_ANNUITIZED-Maintenance_costs) aux_zeros -SpotVector' SpotVector']; %Maximize revenue (Profits-Costs-Investment)
model.objcon=Payment_EFR*Bat_rate*23*365-Transmission_fees*Bat_rate; %23 hours a day available for providing reserve

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
BatterySize2017_EFR(k)=result.x(1);
end
%% Plots
% figure
% plot(Price,BatterySize2017_EFR,'k:','linewidth',0.3);
% hold on 
% plot(Price,BatterySize2017,'k','linewidth',0.3);
% xlim([40 100])
% ylim([0 450]);
% grid
% xlabel('Price (USD/kWh)','Color','k') 
% ylabel('Capacity (MWh)','Color','k')
% legend('With EFR','Without EFR','Interpreter','latex')
% ax = gca;
% ax.FontSize = 12;
% ax.FontName = 'Times';
% set(gcf, 'Color', 'w');
% folder = 'C:\Users\lucas\Documents\Imperial College\Project\Battery investment with EFR\Images';
% baseFileName =  sprintf("EFRComparison2017.pdf");
% fullFileName = fullfile(folder, baseFileName);
% export_fig(fullFileName);
