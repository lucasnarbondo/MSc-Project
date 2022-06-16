% clear all
% clc
% 
% load GurobiDataYear_EFR

%%%%%%%%%%%%%%%%%%%%
%% BATTERIES DATA %%
%%%%%%%%%%%%%%%%%%%%

Bat_investment_cost=90; %USD/kWh
Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000; %USD/MWh/year

Maintenance_costs=0.05*Bat_ANNUITIZED; %Maintenance assumed to be 5% of the capital costs

Bat_rate=50;   %Charge/Discharge rate
Bat_ef=0.94;   %Eficiency

%%%%%%%%%%%%%%%%%%%%%
%% LOAD SPOT PRICE %%
%%%%%%%%%%%%%%%%%%%%%

Year=2018;
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

Transmission_fees=105.38/35*1000*12; %original fees from MIEM
%Transmission_fees=0; %original fees from MIEM

TotalSimulation=length(SpotVector);   %simlulate over a year
aux_zeros=zeros(1,TotalSimulation);   %Auxiliar Vector for algorithm

%%%%%%%%%%%%%
%% PROBLEM %%
%%%%%%%%%%%%%

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
BatterySize=result.x(1);
StorageLevel=result.x(2:1+TotalSimulation);
ChargeRate=result.x(1+TotalSimulation:1+2*TotalSimulation);
DisChargeRate=result.x(1+2*TotalSimulation:1+3*TotalSimulation);

%% Plots
%Plots
% figure
% plot(StorageLevel)
% hold on
% plot(ChargeRate,'--')
% hold on
% plot(DisChargeRate)
% hold on
% plot(SpotVector,':')
% legend('Storage Level','Charge','Discharge','SpotPrice')

fprintf('Optimal Battery Capacity = %.2f \n',result.x(1))