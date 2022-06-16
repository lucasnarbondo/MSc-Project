%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INVESTMENT WITH FREQUENCY RESPONSE %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%%%%%%%%%%%%%%%%%%%%
%% BATTERIES DATA %%
%%%%%%%%%%%%%%%%%%%%

Bat_investment_cost=150; %USD/kWh
Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000/12; %USD/MWh/year

Maintenance_costs=0.05*Bat_ANNUITIZED; %Maintenance assumed to be 5% of the capital costs

Bat_rate=100;   %Charge/Discharge rate
Bat_ef=0.94;   %Eficiency

%%%%%%%%%%%%%%%%%%%%%
%% LOAD SPOT PRICE %%
%%%%%%%%%%%%%%%%%%%%%

Year=2018;
FileName= ['Spot' num2str(Year)];
folder='C:\Users\lucas\Documents\Imperial College\Project\Battery investment with EFR\Spot';
fullFileName = fullfile(folder, FileName);

%SpotPrice=xlsread(fullFileName,'B1:Y365');
SpotPrice=xlsread(fullFileName,'B1:Y31');
SpotVector=reshape(transpose(SpotPrice),[],1); %Re-arrange data as vector


%%%%%%%%%
%% EFR %%
%%%%%%%%%

Time_response=15/60; %15 minutes supplying response
Payment_EFR=12; %Payment for availability to provide EFR in USD/MW/h

%%%%%%%%%%%%%%%%
%% OTHER DATA %%
%%%%%%%%%%%%%%%%

Transmission_fees=105.38/35*1000*12/12; %original fees from MIEM

TotalSimulation=length(SpotVector);   %simlulate over a year
aux_zeros=zeros(1,TotalSimulation);   %Auxiliar Vector for algorithm

%Gurobi Parameters
params.TimeLimit=24*3600;  % Set execution time limit
% Build model
model.modelname  = 'Storage Investment';
model.modelsense = 'max';

%%%%%%%%%%%%%
%% PROBLEM %%
%%%%%%%%%%%%%

%% Define Variables
% index_var=1;
% model.varnames{index_var} = sprintf('Battery Capacity');
% index_var=index_var+1;

% %Battery Storage Level
% for i=1:TotalSimulation
%     model.varnames{index_var} = sprintf('Battery Storage: Time%d',i);
%     index_var=index_var+1;
% end
% 
% %Battery Charge Rate
% for i=1:TotalSimulation
%     model.varnames{index_var} = sprintf('Battery Charge Rate: Time%d',i);
%     index_var=index_var+1;
% end
% 
% %Battery Discharge Rate
% for i=1:TotalSimulation
%     model.varnames{index_var} = sprintf('Battery Discharge Rate: Time%d',i);
%     index_var=index_var+1;
% end


%% Define Objective Function

model.obj=[(-Bat_ANNUITIZED-Maintenance_costs) aux_zeros -SpotVector' SpotVector']; %Maximize revenue (Profits-Costs-Investment)
model.objcon=Payment_EFR*Bat_rate*23*365-Transmission_fees*Bat_rate; %23 hours a day available for providing reserve

%% Build A Matrix - Equality Equations

%Initialize b vector (Ax=b) and sense vector (=,<,>)
b_vector=[];
sense_vector=[];

MatrixIndex=1;

%Storage level Equation: s(t)=s(t-1)+(ef_ch.rate_ch(t)-rate_dis(t)/ef_dis)
for k=1:(TotalSimulation-1)
    aux_A_storage=zeros(1,TotalSimulation);
    aux_A_charge=zeros(1,TotalSimulation);
    aux_A_discharge=zeros(1,TotalSimulation);
    aux_A_storage(k)=1;
    aux_A_storage(k+1)=-1;
    aux_A_charge(k+1)=Bat_ef;
    aux_A_discharge(k+1)=-1/Bat_ef;
    A(MatrixIndex,:)=[0 aux_A_storage aux_A_charge aux_A_discharge];
    model.constrnames{MatrixIndex} = sprintf('Storage Level: Time%d', k+1);
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end

%Initialize Battery at minimum Energy s(t=0)=E_min
E_min=Time_response*Bat_rate;
aux_A_storage=zeros(1,TotalSimulation);
aux_A_storage(1)=1;
A(MatrixIndex,:)=[0 aux_A_storage aux_zeros aux_zeros];
model.constrnames{MatrixIndex} = sprintf('Storage(t=0)=0');
b_vector=[b_vector;E_min];
sense_vector=[sense_vector;'='];
MatrixIndex=MatrixIndex+1;


%% Inequality Equations
%Storage level (t) > E_min     %Needed to have E_min to provide EFR
for k=1:TotalSimulation
    %Create constrains
    aux_A_storage=zeros(1,TotalSimulation);
    aux_A_storage(k)=1;
    b_vector=[b_vector;E_min];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[0 aux_A_storage aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Storage Level > 0: Time%d', k);
    MatrixIndex=MatrixIndex+1;
end

%Storage level (t) < Capacity-E_min
for k=1:TotalSimulation
    %Create constrains
    aux_A_storage=zeros(1,TotalSimulation);
    aux_A_storage(k)=1;
    b_vector=[b_vector;-E_min];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[-1 aux_A_storage aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Storage Level<Capacity-E_min: Time%d', k);
    MatrixIndex=MatrixIndex+1;
end

model.A=sparse(A);
model.rhs=b_vector;
model.sense=sense_vector;


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
figure
plot(StorageLevel)
hold on
plot(ChargeRate,'--')
hold on
plot(DisChargeRate)
hold on
plot(SpotVector,':')
legend('Storage Level','Charge','Discharge','SpotPrice')

fprintf('Optimal Battery Capacity = %.2f \n',result.x(1))