%%%%%%%%%%%%%%%%%%%%%
%% INVESTMENT SPOT %%
%%%%%%%%%%%%%%%%%%%%%

clear all
clc

%%%%%%%%%%%%%%%%%%%%
%% BATTERIES DATA %%
%%%%%%%%%%%%%%%%%%%%

Bat_investment_cost=230; %USD/kWh
Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000/12; %USD/MWh/year (/12 for 1 month test)

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
SpotPrice=xlsread(fullFileName,'B60:Y90');
SpotVector=reshape(transpose(SpotPrice),[],1); %Re-arrange data as vector

%%%%%%%%%%%%%%%%
%% OTHER DATA %%
%%%%%%%%%%%%%%%%

Transmission_fees=105.38/35*1000*12/12; %original fees from MIEM (/12 per month test)

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
model.objcon=-Transmission_fees*Bat_rate; %Transmission fees as a constant

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

%Initialize Battery at 0 s(t=0)=0
aux_A_storage=zeros(1,TotalSimulation);
aux_A_storage(1)=1;
A(MatrixIndex,:)=[0 aux_A_storage aux_zeros aux_zeros];
model.constrnames{MatrixIndex} = sprintf('Storage(t=0)=0');
b_vector=[b_vector;0];
sense_vector=[sense_vector;'='];
MatrixIndex=MatrixIndex+1;


%% Inequality Equations
%Storage level (t) > 0 
for k=1:TotalSimulation
    %Create constrains
    aux_A_storage=zeros(1,TotalSimulation);
    aux_A_storage(k)=1;
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[0 aux_A_storage aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Storage Level > 0: Time%d', k);
    MatrixIndex=MatrixIndex+1;
end

%Storage level (t) < Capacity
for k=1:TotalSimulation
    %Create constrains
    aux_A_storage=zeros(1,TotalSimulation);
    aux_A_storage(k)=1;
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[-1 aux_A_storage aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Storage Level<Capacity: Time%d', k);
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
yyaxis left
plot(StorageLevel2,'linewidth',0.3)
ylim([0 200]);
yyaxis right
plot(ChargeRate,'k:','linewidth',0.3)
hold on
plot(DisChargeRate2,'k--','linewidth',0.3)
yyaxis left
ylabel('StorageLevel (MWh)','Color','k')
yyaxis right
ylabel('Charge/Discharge Rate (MW)','Color','k')
legend('Storage Level','Charge','Discharge','Interpreter','latex')
ylim([0 200]);
grid
xticks([482 494 506 518 530 542 554])
xticklabels({'19-03 00:00','19-03 12:00','20-03 00:00','20-03 12:00','21-03 00:00','21-03 12:00','22-03 00:00'})
xtickangle(45)
pbaspect([2 1 1])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');

figure
yyaxis left
plot(StorageLevel,'linewidth',0.3)
ylim([0 200]);
hold on
yyaxis right
plot(SpotVector,'k:','linewidth',0.3)
yyaxis left
ylabel('Storage Level (MWh)','Color','k')
yyaxis right
ylabel('Spot Price (USD/MWh)','Color','k')
legend('Storage Level','Spot Price','Interpreter','latex')
grid
xticks([482 494 506 518 530 542 554])
xticklabels({'19-03 00:00','19-03 12:00','20-03 00:00','20-03 12:00','21-03 00:00','21-03 12:00','22-03 00:00'})
xtickangle(45)
pbaspect([2 1 1])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');

fprintf('Optimal Battery Capacity = %.2f \n',result.x(1))