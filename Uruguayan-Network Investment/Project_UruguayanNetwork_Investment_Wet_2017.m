%%%%%%%%%%%%%%%%%%
%% THREE BUS TEST%
%%%%%%%%%%%%%%%%%%

clear all
clc

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
load GenerationAndDemandData_Wet_2017
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

Bat_investment_cost=200; %USD/kWh
Interest_rate=0.05; % 5% interest rate
number_years=15;
Ann_factor=(1-(1/(1+Interest_rate)^number_years))/Interest_rate; %Factor for annuitized value

Bat_ANNUITIZED=Bat_investment_cost/Ann_factor*1000; %USD/MWh/year
Maintenance_costs=0.05*Bat_ANNUITIZED; %Maintenance assumed to be 5% of the capital costs

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

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FORMULATE THE PROBLEM %
%%%%%%%%%%%%%%%%%%%%%%%%%%

HourlyCost=[];
LIIONStorageVar=[];
RateChargeLIIONVar=[];
RateDisChargeLIIONVar=[];
AngleDeltaVar=[];
StorageNodeK=zeros(1,NumberBuses);


%% Define Variables
%Production
index_var=1;
% for k=1:TotalGenerators
%     for i=1:TotalSimulation
%         model.varnames{index_var} = sprintf('Thermal: Number%d Time%d',k,i);
%         index_var=index_var+1;
%     end
% end
% 
% %LI-ION Storage Leve
% for k=1:NumberBuses
%     for i=1:TotalSimulation
%         model.varnames{index_var} = sprintf('LI-ION Storage: Bus%d Time%d',k,i);
%         index_var=index_var+1;
%     end
% end
% 
% %LI-ION Charge Rate
% for k=1:NumberBuses
%     for i=1:TotalSimulation
%         model.varnames{index_var} = sprintf('Li-ION Charge Rate: Bus%d Time%d',k,i);
%         index_var=index_var+1;
%     end
% end
% 
% %LI-ION Discharge Rate
% for k=1:NumberBuses
%     for i=1:TotalSimulation
%         model.varnames{index_var} = sprintf('LI-ION Discharge Rate: Bus%d Time%d',k,i);
%         index_var=index_var+1;
%     end
% end
% 
% %Angle delta
% for k=1:NumberBuses
%     for i=1:TotalSimulation
%         model.varnames{index_var} = sprintf('Angle: Bus%d Time%d',k,i);
%         index_var=index_var+1;
%     end
% end
% 
% %K-Siting LI-ION
% for k=1:NumberBuses
%         model.varnames{index_var} = sprintf('K siting LI-ION: Bus%d',k);
%         index_var=index_var+1;
% end

%% Objective Function
aux_zeros=zeros(1,NumberBuses*TotalSimulation);
aux_zeros_gen=zeros(1,2*TotalSimulation);

GenerationProduction=ones(1,TotalSimulation);    %Production vector 
HourlyCost=[FuelCost1*GenerationProduction FuelCost2*GenerationProduction];     %Hourly costs for objective function

StorageNodeK_invest=(Bat_ANNUITIZED+Maintenance_costs)*ones(1,NumberBuses)*100; %*100 because of per unit

%Define Objective Function: Minimize Production Cost
%Rest of the variables for optimization constrains
model.obj=[HourlyCost*52 aux_zeros aux_zeros aux_zeros StorageNodeK_invest aux_zeros aux_zeros aux_zeros];


%% Build A Matrix - Equality Equations

%Initialize b vector (Ax=b) and sense vector (=,<,>)
b_vector=[];
sense_vector=[];
MatrixIndex=1;

%Storage level Equation: s(t)=s(t-1)+(ef_ch.rate_ch(t)-rate_dis(t)/ef_dis)deltaT
for n=1:NumberBuses
for k=1:(TotalSimulation-1)
    %Create constrains for Li_ION Storage
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_charge_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    aux_A_storage_LIION(TotalSimulation*(n-1)+k+1)=-1;
    aux_A_charge_LIION(TotalSimulation*(n-1)+k+1)=LIION_ef*TimeStep;
    A(MatrixIndex,:)=[aux_zeros_gen aux_A_storage_LIION aux_A_charge_LIION aux_zeros StorageNodeK aux_zeros aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('LI-ION Storage Level: Bus%d Time%d', n, k+1);
    %Create constrains for batteries storage
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end
end

%Initialize Battery at 0 s(t=0)=0
for n=1:NumberBuses
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+1)=1;
    A(MatrixIndex,:)=[aux_zeros_gen aux_A_storage_LIION aux_zeros aux_zeros StorageNodeK aux_zeros aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('LiION Storage(0): Bus%d', n);
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end


%Cyclic storage level Equation: s(t=1)=s(t=T)
%Battery
% for n=1:NumberBuses
%     aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
%     aux_A_storage_LIION(TotalSimulation*(n-1)+1)=1;
%     aux_A_storage_LIION(TotalSimulation*(n-1)+TotalSimulation)=-1;
%     A(MatrixIndex,:)=[aux_zeros_gen aux_A_storage_LIION aux_zeros aux_zeros StorageNodeK aux_zeros aux_zeros aux_zeros];
%     model.constrnames{MatrixIndex} = sprintf('Cyclic LiION Storage: Bus%d', n);
%     b_vector=[b_vector;0];
%     sense_vector=[sense_vector;'='];
%     MatrixIndex=MatrixIndex+1;
% end

%Construct Adjacent Matrix for nodes
ADJ=sparse(Ob', Db', ones((length(Db)),1), length(Ob+1), length(Ob+1));
ADJ_X=sparse(Ob', Db', xr, length(Ob+1), length(Ob+1));
ADJ_X = (ADJ_X+ADJ_X') - eye(size(ADJ_X,1)).*diag(ADJ_X); %Make ADJ_X Symmetrical
Demand=0;
Wind=0;
ExpArg=0;
ExpBrasil=0;
Solar=0;
Biomass=0;

%Power Balance Equations:
for i=1:NumberBuses     %Number of equations: TotalSimulation*NumberBuses
    for m=1:TotalSimulation
        for aux=1:NumberBuses %Reset Angle Delta for each equation
        AngleDelta{aux}=zeros(1,TotalSimulation); 
        end
        
        for j=1:NumberBuses %Verification if node j is connected to node i
            if (ADJ(i,j)==1)||(ADJ(j,i)==1)  %If nodes are connected:
                AngleDelta{i}(m)=AngleDelta{i}(m)-(-1/ADJ_X(i,j));      %Set value for equation (Delta_n)
                AngleDelta{j}(m)=-1/ADJ_X(i,j);                         %Set value for equation (Delta_m)
            end
        
        aux_A_delta=[];    
        for aux=1:NumberBuses
        aux_A_delta=[aux_A_delta AngleDelta{aux}];
        end
            
        %RATE CHARGE-DISCHARGE LI-ION
        aux_A_charge_LIION=zeros(1,NumberBuses*TotalSimulation);
        aux_A_charge_LIION(TotalSimulation*(i-1)+m)=-1;
        
        %WIND CURTAILMENT
        aux_A_curtailmentWind=zeros(1,NumberBuses*TotalSimulation);
        aux_A_curtailmentWind(TotalSimulation*(i-1)+m)=1;
        
        %HYDRO CURTAILMENT
        aux_A_curtailmentHydroSG=zeros(1,NumberBuses*TotalSimulation);
        aux_A_curtailmentHydroSG(TotalSimulation*(i-1)+m)=1;
        aux_A_curtailmentHydroRN=zeros(1,NumberBuses*TotalSimulation);
        aux_A_curtailmentHydroRN(TotalSimulation*(i-1)+m)=1;        
        
        %PRODUCTION NODE i
        aux_A_production_1=zeros(1,TotalSimulation);
        aux_A_production_2=zeros(1,TotalSimulation);
        if Generator1Bus(i)==1
            aux_A_production_1(m)=1;
        end
        if Generator2Bus(i)==1
            aux_A_production_2(m)=1;
        end
                
        %DEMAND (t)
        Demand=DemandAvg(m)*DemandBus(i);
        %EXP ARGENTINA (t)
        ExpArg=ExpArgAvg(m)*ExpArgBus(i);
        %EXP BRASIL (t)
        ExpBrasil=ExpBrasilAvg(m)*ExpBrasilBus(i);
        
        %WIND(t)
        Wind=WindAvg(m)*WindBus(i);
        %SOLAR(t)
        Solar=SolarAvg(m)*SolarBus(i);
        %BIOMASS(t)
        Biomass=BiomassAvg(m)*BiomassBus(i);
        %SaltoGrande(t)
        SaltoGrande=SaltoGrandeAvg(m)*SaltoGrandeBus(i);
        %RioNegro(t)
        RioNegro=RioNegroAvg(m)*RioNegroBus(i);
        end
        
    A(MatrixIndex,:)=[aux_A_production_1 aux_A_production_2 aux_zeros aux_A_charge_LIION aux_A_delta StorageNodeK aux_A_curtailmentWind aux_A_curtailmentHydroSG aux_A_curtailmentHydroRN];
    model.constrnames{MatrixIndex} = sprintf('Power Balance: Bus%d Time%d', i, m);
    MatrixIndex=MatrixIndex+1;  
    b_vector=[b_vector;Demand+ExpArg+ExpBrasil-Wind-Solar-Biomass-SaltoGrande-RioNegro];
    sense_vector=[sense_vector;'='];
    end
end

%Finally set slack bus angle = 0
for i=1:TotalSimulation
    SwingBus=1;
    aux_A_delta=zeros(1,NumberBuses*TotalSimulation);
    aux_A_delta(TotalSimulation*(SwingBus-1)+i)=1;
    A(MatrixIndex,:)=[aux_zeros_gen aux_zeros aux_zeros aux_A_delta StorageNodeK aux_zeros aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Swing Angle: Time%d', i);
    MatrixIndex=MatrixIndex+1;  
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
end


%% A Matrix - Inequality Equations

%Storage level (t) > 0 LI-ION
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for LI-ION Storage
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[aux_zeros_gen aux_A_storage_LIION aux_zeros aux_zeros StorageNodeK aux_zeros aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Storage Level LI-ION >0: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end


%Storage level (t) < Kj LI-ION
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for LI-ION
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    StorageNodeK=zeros(1,NumberBuses);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    StorageNodeK(n)=-1; 
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[aux_zeros_gen aux_A_storage_LIION aux_zeros aux_zeros StorageNodeK aux_zeros aux_zeros aux_zeros];
    model.constrnames{MatrixIndex} = sprintf('Storage Level LI-ION <kj: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end

StorageNodeK=zeros(1,NumberBuses);
%Charge<WindCurtailed+HydroCurtailed
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for LI-ION
    aux_A_charge_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_WindCurtailed=zeros(1,NumberBuses*TotalSimulation);
    aux_A_HydroSGCurtailed=zeros(1,NumberBuses*TotalSimulation);
    aux_A_HydroRNCurtailed=zeros(1,NumberBuses*TotalSimulation);
    aux_A_charge_LIION(TotalSimulation*(n-1)+k)=1;
    aux_A_WindCurtailed(TotalSimulation*(n-1)+k)=-1; 
    aux_A_HydroSGCurtailed(TotalSimulation*(n-1)+k)=-1; 
    aux_A_HydroRNCurtailed(TotalSimulation*(n-1)+k)=-1; 
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[aux_zeros_gen aux_zeros aux_A_charge_LIION aux_zeros StorageNodeK aux_A_WindCurtailed aux_A_HydroSGCurtailed aux_A_HydroRNCurtailed];
    model.constrnames{MatrixIndex} = sprintf('Storage Level LI-ION <kj: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end


% %Ramp Rate Generators > -RR
% for n=1:TotalGenerators
% for k=1:TotalSimulation-1
%     aux_A_production=zeros(1,TotalGenerators*TotalSimulation);
%     aux_A_production(TotalSimulation*(n-1)+k)=-1;
%     aux_A_production(TotalSimulation*(n-1)+k+1)=1;
%     b_vector=[b_vector;-RampRate];
%     sense_vector=[sense_vector;'>'];
%     A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
%     model.constrnames{MatrixIndex} = sprintf('p(t)-p(t-1)>-RR: Bus%d Time%d', n, k);
%     MatrixIndex=MatrixIndex+1;
% end
% end
% %Ramp Rate Generators < RR
% for n=1:TotalGenerators
% for k=1:TotalSimulation-1
%     aux_A_production=zeros(1,TotalGenerators*TotalSimulation);
%     aux_A_production(TotalSimulation*(n-1)+k)=-1;
%     aux_A_production(TotalSimulation*(n-1)+k+1)=1;
%     b_vector=[b_vector;RampRate];
%     sense_vector=[sense_vector;'<'];
%     A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
%     model.constrnames{MatrixIndex} = sprintf('p(t)-p(t-1)<RR: Bus%d Time%d', n, k);
%     MatrixIndex=MatrixIndex+1;
% end
% end
% 

%Flow through line > -TCmax
ADJ_capacity=sparse(Ob', Db', F0, length(Ob+1), length(Ob+1));
ADJ_capacity = (ADJ_capacity+ADJ_capacity') - eye(size(ADJ_capacity,1)).*diag(ADJ_capacity); %Make ADJ_capacity Symmetrical
ADJ_capacity=full(ADJ_capacity);
for m=1:TotalSimulation
for i=1:NumberBuses
     for j=1:NumberBuses %Verification if node j is connected to node i
        for aux=1:NumberBuses %Reset Angle Delta for each equation
            AngleDelta{aux}=zeros(1,TotalSimulation); 
        end
            if (ADJ(i,j)==1)||(ADJ(j,i)==1)  %If nodes are connected:
                AngleDelta{i}(m)=(-1/ADJ_X(i,j));      %Set value for equation (Delta_n)
                AngleDelta{j}(m)=1/ADJ_X(i,j);         %Set value for equation (Delta_m)
                TCmax=ADJ_capacity(i,j);
                aux_A_delta=[];    
                    for aux=1:NumberBuses
                    aux_A_delta=[aux_A_delta AngleDelta{aux}];
                    end
                b_vector=[b_vector;-TCmax];
                sense_vector=[sense_vector;'>'];
                A(MatrixIndex,:)=[aux_zeros_gen aux_zeros aux_zeros aux_A_delta StorageNodeK aux_zeros aux_zeros aux_zeros];
                model.constrnames{MatrixIndex} = sprintf('Flow>-TCmax: Line(%d,%d) Time%d', i, j, m);
                MatrixIndex=MatrixIndex+1;
            end

     end

end
end 

% %Flow Through Line < TCmax
for m=1:TotalSimulation
for i=1:NumberBuses
     for j=1:NumberBuses %Verification if node j is connected to node i
        for aux=1:NumberBuses %Reset Angle Delta for each equation
            AngleDelta{aux}=zeros(1,TotalSimulation); 
        end
            if (ADJ(i,j)==1)||(ADJ(j,i)==1)  %If nodes are connected:
                AngleDelta{i}(m)=(-1/ADJ_X(i,j));      %Set value for equation (Delta_n)
                AngleDelta{j}(m)=1/ADJ_X(i,j);         %Set value for equation (Delta_m)
                TCmax=ADJ_capacity(i,j);
                aux_A_delta=[];    
                    for aux=1:NumberBuses
                    aux_A_delta=[aux_A_delta AngleDelta{aux}];
                    end
                b_vector=[b_vector;TCmax];
                sense_vector=[sense_vector;'<'];
                A(MatrixIndex,:)=[aux_zeros_gen aux_zeros aux_zeros aux_A_delta StorageNodeK aux_zeros aux_zeros aux_zeros];
                model.constrnames{MatrixIndex} = sprintf('Flow<TCmax: Line(%d,%d) Time%d', i, j, m);
                MatrixIndex=MatrixIndex+1;
            end

     end

end
end

model.A=sparse(A);
model.rhs=b_vector;
model.sense=sense_vector;


%% Bounded Equations

LowerBound=[];
UpperBound=[];

%Production
 for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;GeneratorCapacity1];
 end
 for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;GeneratorCapacity2];
 end
    

%LI-ION Storage Leve
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;inf];
    end
end
%LI-ION Charge Rate
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;-LIION_DR];
    UpperBound=[UpperBound;LIION_DR];
    end
end
%LI-ION Discharge Rate

%Angle delta
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;-pi];
    UpperBound=[UpperBound;pi];
    end
end
%K-Siting LI-ION
for i=1:NumberBuses
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;inf];
end

%Wind Curtailed
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    if ThermalAvg(k) > 0
        UpperBound=[UpperBound;0];
    else
        UpperBound=[UpperBound;WindBus(i)*4];
    end
    end
end

%HydroSG Curtailed
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    if ThermalAvg(k) > 0
        UpperBound=[UpperBound;0];
    else
        UpperBound=[UpperBound;(9.45-SaltoGrandeAvg(k))*SaltoGrandeBus(i)];
    end
    end
end

%HydroRN Curtailed
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    if ThermalAvg(k) > 0
        UpperBound=[UpperBound;0];
    elseif 3.33+1.6-RioNegroAvg(k)>0
        UpperBound=[UpperBound;(3.33+1.6-RioNegroAvg(k))*RioNegroBus(i)];
        else
        UpperBound=[UpperBound;0];
    end
    end
end

model.lb=LowerBound;
model.ub=UpperBound;


%% Gurobi Optimization
tic
result=gurobi(model,params);
toc

%% Plot Results
ProductionThermal1=result.x(1:TotalSimulation);
ProductionThermal2=result.x(TotalSimulation+1:2*TotalSimulation);
LiIONStorageLevelBus1=result.x(2*TotalSimulation+1:3*TotalSimulation);
LiIONStorageLevelBus2=result.x(3*TotalSimulation+1:4*TotalSimulation);
LiIONStorageLevelBus3=result.x(4*TotalSimulation+1:5*TotalSimulation);
LiIONStorageLevelBus4=result.x(5*TotalSimulation+1:6*TotalSimulation);
LiIONStorageLevelBus5=result.x(6*TotalSimulation+1:7*TotalSimulation);
LiIONStorageLevelBus6=result.x(7*TotalSimulation+1:8*TotalSimulation);
LiIONStorageLevelBus7=result.x(8*TotalSimulation+1:9*TotalSimulation);
LiIONChargeRateBus1 = result.x(9*TotalSimulation+1:10*TotalSimulation);
LiIONChargeRateBus2 = result.x(10*TotalSimulation+1:11*TotalSimulation);
LiIONChargeRateBus3 = result.x(11*TotalSimulation+1:12*TotalSimulation);
LiIONChargeRateBus4 = result.x(12*TotalSimulation+1:13*TotalSimulation);
LiIONChargeRateBus5 = result.x(13*TotalSimulation+1:14*TotalSimulation);
LiIONChargeRateBus6 = result.x(14*TotalSimulation+1:15*TotalSimulation);
LiIONChargeRateBus7 = result.x(15*TotalSimulation+1:16*TotalSimulation);

fprintf('Li-ION Siting Bus1 = %.2f \n',result.x(23*TotalSimulation+1)*100)
fprintf('Li-ION Siting Bus2 = %.2f \n',result.x(23*TotalSimulation+2)*100)
fprintf('Li-ION Siting Bus3 = %.2f \n',result.x(23*TotalSimulation+3)*100)
fprintf('Li-ION Siting Bus4 = %.2f \n',result.x(23*TotalSimulation+4)*100)
fprintf('Li-ION Siting Bus5 = %.2f \n',result.x(23*TotalSimulation+5)*100)
fprintf('Li-ION Siting Bus6 = %.2f \n',result.x(23*TotalSimulation+6)*100)
fprintf('Li-ION Siting Bus7 = %.2f \n',result.x(23*TotalSimulation+7)*100)

WindCurtailmentBus1 = result.x(23*TotalSimulation+8:24*TotalSimulation+7);
WindCurtailmentBus2 = result.x(24*TotalSimulation+8:25*TotalSimulation+7);
WindCurtailmentBus3 = result.x(25*TotalSimulation+8:26*TotalSimulation+7);
WindCurtailmentBus4 = result.x(26*TotalSimulation+8:27*TotalSimulation+7);
WindCurtailmentBus5 = result.x(27*TotalSimulation+8:28*TotalSimulation+7);
WindCurtailmentBus6 = result.x(28*TotalSimulation+8:29*TotalSimulation+7);
WindCurtailmentBus7 = result.x(29*TotalSimulation+8:30*TotalSimulation+7);

SGCurtailment = result.x(30*TotalSimulation+8:31*TotalSimulation+7);

RNCurtailment = result.x(39*TotalSimulation+8:40*TotalSimulation+7);

DemandAvg=DemandAvg(1:TotalSimulation);
ExpArgAvg=ExpArgAvg(1:TotalSimulation);
ExpBrasilAvg=ExpBrasilAvg(1:TotalSimulation);
WindAvg=WindAvg(1:TotalSimulation);
SolarAvg=SolarAvg(1:TotalSimulation);
BiomassAvg=BiomassAvg(1:TotalSimulation);
ThermalAvg=ThermalAvg(1:TotalSimulation);
SaltoGrandeAvg=SaltoGrandeAvg(1:TotalSimulation);
RioNegroAvg=RioNegroAvg(1:TotalSimulation);

ThermalProduction=ProductionThermal1+ProductionThermal2;
TotalStorageLevel=LiIONStorageLevelBus1+LiIONStorageLevelBus2+LiIONStorageLevelBus3+LiIONStorageLevelBus4+LiIONStorageLevelBus5+LiIONStorageLevelBus6+LiIONStorageLevelBus7;
TotalChargeRate=LiIONChargeRateBus1+LiIONChargeRateBus2+LiIONChargeRateBus3+LiIONChargeRateBus4+LiIONChargeRateBus5+LiIONChargeRateBus6+LiIONChargeRateBus7;
TotalWindCurtailed=WindCurtailmentBus1+WindCurtailmentBus2+WindCurtailmentBus3+WindCurtailmentBus4+WindCurtailmentBus5+WindCurtailmentBus6+WindCurtailmentBus7;


%Plot Overall System
figure
plot((DemandAvg+ExpArgAvg+ExpBrasilAvg)*100,'k','linewidth',0.3)
hold on
plot((WindAvg+SolarAvg+BiomassAvg)*100,'k--','linewidth',0.3)
hold on
plot((SaltoGrandeAvg+RioNegroAvg)*100,'k:','linewidth',0.3)
hold on
plot(ThermalProduction*100,'linewidth',0.3)
hold on
plot(TotalStorageLevel*100,'linewidth',0.3)
hold on
plot(TotalWindCurtailed*100,'linewidth',0.3)
hold on
plot(SGCurtailment*100,'linewidth',0.3)
hold on
plot(RNCurtailment*100,'linewidth',0.3)
legend('Demand+Exportation','Wind+Solar+Biomass','Hydro','Thermal Production','Storage Level','Wind Previously Curtailed','SG Previously Curtailed','RN Previously Curtailed','Location','northwest','Interpreter','latex')
grid
xticks([12 36 60 84 108 132 156])
xticklabels({'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7'})
xtickangle(45)
ylabel('Energy (MWh)','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network Investment\Images';
baseFileName =  sprintf("All2017Wet-Case1.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);

%Plot Storage

figure
yyaxis left
plot(TotalStorageLevel*100,'linewidth',0.3)
yyaxis right
plot(TotalChargeRate*100,'--','linewidth',0.3)
yyaxis left
ylabel('StorageLevel (MWh)','Color','k')
yyaxis right
ylabel('Charge/Discharge Rate (MW)','Color','k')
legend('Storage Level','Charge','Interpreter','latex')
grid
xticks([12 36 60 84 108 132 156])
xticklabels({'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7'})
xtickangle(45)
pbaspect([2 1 1])
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network Investment\Images';
baseFileName =  sprintf("ChargeandDischargeWet2017-Case1.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);

%Plot Thermal
figure
plot(ThermalAvg*100,'k','linewidth',0.3)
hold on
plot(ThermalProduction*100,'k:','linewidth',0.3)
legend('Thermal pre-battery','Thermal post-battery','Location','northwest','Interpreter','latex')
xticks([12 36 60 84 108 132 156])
xticklabels({'Day 1','Day 2','Day 3','Day 4','Day 5','Day 6','Day 7'})
xtickangle(45)
ylabel('Energy (MWh)','Color','k')
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Times';
set(gcf, 'Color', 'w');
folder = 'C:\Users\lucas\Documents\Imperial College\Project\Uruguayan Network Investment\Images';
baseFileName =  sprintf("ThermalPrePostWet2017-Case1.pdf");
fullFileName = fullfile(folder, baseFileName);
export_fig(fullFileName);

%Plot Error
% figure
% TotalDemand=DemandAvg+ExpArgAvg+ExpBrasilAvg+TotalChargeRate;
% TotalPower=WindAvg+SolarAvg+BiomassAvg+ThermalProduction+SaltoGrandeAvg+RioNegroAvg+TotalWindCurtailed+SGCurtailment+RNCurtailment;
% plot(TotalDemand-TotalPower)




% PDT=zeros(1,TotalSimulation);
% CTR=zeros(1,TotalSimulation);
% 
% for k=1:TotalSimulation
%     if ThermalAvg(k)-3.6>0
%         CTR(k)=ThermalAvg(k)-3.6;
%         PDT(k)=3.6;
%     else
%         CTR(k)=0;
%         PDT(k)=ThermalAvg(k);
%     end
% end

