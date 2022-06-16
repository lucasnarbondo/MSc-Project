%%%%%%%%%%%%%%%%%%
%% THREE BUS TEST%
%%%%%%%%%%%%%%%%%%

clear all
clc
load GurobiDataWeek2017_Dry_Case2


%%%%%%%%%%%%%%%%
% NETWORK DATA %
%%%%%%%%%%%%%%%%

NumberBuses=7;  % number of buses 1)SG 2)San Javier 3)Palmar 4)MontevideoA 5)MontevideoB 6)San Carlos 7)Melo

xr=[0.05 0.05 0.05 0.05 0.05 0.05 0.05]; % reactance of the line [p.u.]
Ob=[1 2 3 3 4 4 6];       % Origin bus
Db=[2 3 4 5 5 6 7];       % Destination bus
F0=[1380, 1380, 1380, 1380, 1380, 1380, 1380]./100;    % Lines capacity (MW)
%F0=[999999, 999999, 99999, 999999, 99999, 9999999, 9999999]./100;    % Lines capacity (MW)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THERMAL GENERATORS DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

TotalGenerators=2;
GeneratorCapacity1=280/100;   % Thermal generator's capacity (MW) 
FuelCost1=140*100;            % Operation Cost £/MWh
Generator1Bus=[0 0 0 1 0 0 0];
GeneratorCapacity2=360/100;
FuelCost2=120*100;
Generator2Bus=[0 0 0 0 1 0 0];
ThermalGeneratorBus=[4 5];

%%%%%%%%%%%%%%%%%%%%%%%
% WIND AND SOLAR DATA %
%%%%%%%%%%%%%%%%%%%%%%%

%Wind data over a day
load GenerationAndDemandData_Dry_2017
WindAvg=WindAvg/100; %p.u
SolarAvg=SolarAvg/100;


%Wind in each bus
WindBus=[0.10 0 0.20 0.20 0.20 0.2 0.10];   
SolarBus=[1 0 0 0 0 0 0];

%%%%%%%%%%%%%%%%%%%
% HIDROPOWER DATA %
%%%%%%%%%%%%%%%%%%%
SaltoGrandeAvg=SaltoGrandeAvg/100;
SaltoGrandeBus=[1 0 0 0 0 0 0];

RioNegroAvg=RioNegroAvg/100;
RioNegroBus=[0 0 1 0 0 0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BIOMASS AND THERMAL DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiomassAvg=BiomassAvg/100;
BiomassBus=[0 1 0 0 0 0 0];

ThermalAvg=ThermalAvg/100;

%%%%%%%%%%%%%%%%
% STORAGE DATA %
%%%%%%%%%%%%%%%%

Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

LIION_ef=0.94;
LIION_DR=100/100;

%%%%%%%%%%%%%%
% OTHER DATA %
%%%%%%%%%%%%%%

TimeStep=60/60;                        %hours
TotalSimulation=24*7/TimeStep;         %simlulate over a week 
RampRate=4*6/100;                      %MW/time step

%Demand data
DemandBus=[0.05 0.05 0.05 0.30 0.30 0.20 0.05];
DemandAvg=DemandAvg/100; %p.u

ExpArgBus=[1 0 0 0 0 0 0];
ExpArgAvg=ExpArgAvg/100; %p.u

ExpBrasilBus=[0 0 0 0 0 0 1];
ExpBrasilAvg=ExpBrasilAvg/100; %p.u

%Gurobi Parameters
params.TimeLimit=24*3600;  % Set execution time limit
% Build model
model.modelname  = 'Storage Siting';
model.modelsense = 'min';

%%%%%%%%%%%%%
%% PROBLEM %%
%%%%%%%%%%%%%

BatterySize=zeros(1,100);
Price=zeros(1,100);

for k=1:100
Bat_investment_cost=40+1.5*(k-1);
Price(k)=Bat_investment_cost;

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000; %USD/MWh/year

Maintenance_costs=0.05*Bat_ANNUITIZED;


%% Objective Function
aux_zeros=zeros(1,NumberBuses*TotalSimulation);
aux_zeros_gen=zeros(1,2*TotalSimulation);

GenerationProduction=ones(1,TotalSimulation);    %Production vector 
HourlyCost=[FuelCost1*GenerationProduction FuelCost2*GenerationProduction];     %Hourly costs for objective function

StorageNodeK_invest=(Bat_ANNUITIZED+Maintenance_costs)*ones(1,NumberBuses)*100; %*100 because of per unit

%Define Objective Function: Minimize Production Cost
%Rest of the variables for optimization constrains
model.obj=[HourlyCost*52 aux_zeros aux_zeros aux_zeros StorageNodeK_invest aux_zeros];

%% Gurobi Optimization
tic
result=gurobi(model,params);
toc

BatteryBus1=result.x(23*TotalSimulation+1)*100;
BatteryBus2=result.x(23*TotalSimulation+2)*100;
BatteryBus3=result.x(23*TotalSimulation+3)*100;
BatteryBus4=result.x(23*TotalSimulation+4)*100;
BatteryBus5=result.x(23*TotalSimulation+5)*100;
BatteryBus6=result.x(23*TotalSimulation+6)*100;
BatteryBus7=result.x(23*TotalSimulation+7)*100;
TotalBatterySize=BatteryBus1+BatteryBus2+BatteryBus3+BatteryBus4+BatteryBus5+BatteryBus6+BatteryBus7;

BatterySize(k)=TotalBatterySize;
end

figure
plot(Price,BatterySize,'k','linewidth',0.3);
grid
xlabel('Price (USD/kWh)','Color','k') 
ylabel('Capacity (MWh)','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network Investment\Images';
baseFileName =  sprintf("PriceVariationDry2017-Case2.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);
