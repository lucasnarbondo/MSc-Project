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

%%%%%%%%%%%%%%%%%%
% GENERATOR DATA %
%%%%%%%%%%%%%%%%%%

NumberGenerators=3;                           % number of installed generators
%GeneratorCapacity=[945 0 593 280 350 0 0]./100;   % generator's capacity (MW) %SG and Rio Negro installed capacity
GeneratorCapacity=[450 10 360 280 350 10 10]./100;
FuelCost=[0 0 0 140 120 0 0].*100;                     % Operation Cost £/MWh


%%%%%%%%%%%%%%%%%%%%%%%
% WIND AND SOLAR DATA %
%%%%%%%%%%%%%%%%%%%%%%%

%Wind data over a day
load GenerationAndDemandData2
WindAvg=WindAvg/100; %p.u
SolarAvg=SolarAvg/100;


%Wind in each bus
WindBus=[0.10 0 0.20 0.20 0.20 0.2 0.10];   
SolarBus=[1 0 0 0 0 0 0];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BIOMASS AND THERMAL DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiomassAvg=BiomassAvg/100;
BiomassBus=[0 1 0 0 0 0 0];

ThermalAvg=ThermalAvg/100;

%%%%%%%%%%%%%%%%
% STORAGE DATA %
%%%%%%%%%%%%%%%%

%LI-ION Batteries
LIION_ef=0.94;
LIION_SC=100/100;
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
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('Production: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%LI-ION Storage Leve
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('LI-ION Storage: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%LI-ION Charge Rate
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('Li-ION Charge Rate: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%LI-ION Discharge Rate
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('LI-ION Discharge Rate: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%Angle delta
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('Angle: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%K-Siting LI-ION
for k=1:NumberBuses
        model.varnames{index_var} = sprintf('K siting LI-ION: Bus%d',k);
        index_var=index_var+1;
end

%% Objective Function
for k=1:NumberBuses
    GenerationProduction{k}=ones(1,TotalSimulation);    %Production vector 
    HourlyCost=[HourlyCost FuelCost(k)*GenerationProduction{k}];     %Hourly costs for objective function
    
    %LI-ION Batteries
    LIIONStorage{k}=zeros(1,TotalSimulation);
    LIIONStorageVar=[LIIONStorageVar 0*LIIONStorage{k}];                        %Variable for optimization problem
    
    %Rate Charge LIION
    RateChargeLIION{k}=zeros(1,TotalSimulation); 
    RateChargeLIIONVar=[RateChargeLIIONVar 0*RateChargeLIION{k}];           %Variable for optimization problem 
    
    %Rate DisCharge LIION
    RateDisChargeLIION{k}=zeros(1,TotalSimulation); 
    RateDisChargeLIIONVar=[RateDisChargeLIIONVar 0*RateDisChargeLIION{k}];  %Variable for optimization problem

    %Angles
    AngleDelta{k}=zeros(1,TotalSimulation); 
    AngleDeltaVar=[AngleDeltaVar 0*AngleDelta{k}];  %Variable for optimization problem
end

%Define Objective Function: Minimize Production Cost
%Rest of the variables for optimization constrains
model.obj=[HourlyCost LIIONStorageVar RateChargeLIIONVar RateDisChargeLIIONVar AngleDeltaVar StorageNodeK];


%% Build A Matrix - Equality Equations

aux_zeros=zeros(1,NumberBuses*TotalSimulation);

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
    aux_A_discharge_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    aux_A_storage_LIION(TotalSimulation*(n-1)+k+1)=-1;
    aux_A_charge_LIION(TotalSimulation*(n-1)+k+1)=LIION_ef*TimeStep;
    aux_A_discharge_LIION(TotalSimulation*(n-1)+k+1)=-TimeStep/LIION_ef;
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_LIION aux_A_charge_LIION aux_A_discharge_LIION aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('LI-ION Storage Level: Bus%d Time%d', n, k+1);
    %Create constrains for batteries storage
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end
end

%Cyclic storage level Equation: s(t=1)=s(t=T)
%Battery
for n=1:NumberBuses
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+1)=1;
    aux_A_storage_LIION(TotalSimulation*(n-1)+TotalSimulation)=-1;
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Cyclic LiION Storage: Bus%d', n);
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end

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
        aux_A_discharge_LIION=zeros(1,NumberBuses*TotalSimulation);
        aux_A_charge_LIION(TotalSimulation*(i-1)+m)=-1;
        aux_A_discharge_LIION(TotalSimulation*(i-1)+m)=1;
        
        %PRODUCTION NODE i
        aux_A_production=zeros(1,NumberBuses*TotalSimulation);
        aux_A_production(TotalSimulation*(i-1)+m)=1;
        
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
        end
        
    A(MatrixIndex,:)=[aux_A_production aux_zeros aux_A_charge_LIION aux_A_discharge_LIION aux_A_delta StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Power Balance: Bus%d Time%d', i, m);
    MatrixIndex=MatrixIndex+1;  
    b_vector=[b_vector;Demand+ExpArg+ExpBrasil-Wind-Solar-Biomass];
    sense_vector=[sense_vector;'='];
    end
end

%Finally set slack bus angle = 0
for i=1:TotalSimulation
    [maxCapacity, SwingBus]= max(GeneratorCapacity);
    aux_A_delta=zeros(1,NumberBuses*TotalSimulation);
    aux_A_delta(TotalSimulation*(SwingBus-1)+i)=1;
    A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_A_delta StorageNodeK];
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
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_zeros StorageNodeK];
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
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Storage Level LI-ION <kj: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end
StorageNodeK=zeros(1,NumberBuses);

%Ramp Rate Generators > -RR
for n=1:NumberBuses
for k=1:TotalSimulation-1
    aux_A_production=zeros(1,NumberBuses*TotalSimulation);
    aux_A_production(TotalSimulation*(n-1)+k)=-1;
    aux_A_production(TotalSimulation*(n-1)+k+1)=1;
    b_vector=[b_vector;-RampRate];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('p(t)-p(t-1)>-RR: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end
%Ramp Rate Generators < RR
for n=1:NumberBuses
for k=1:TotalSimulation-1
    aux_A_production=zeros(1,NumberBuses*TotalSimulation);
    aux_A_production(TotalSimulation*(n-1)+k)=-1;
    aux_A_production(TotalSimulation*(n-1)+k+1)=1;
    b_vector=[b_vector;RampRate];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('p(t)-p(t-1)<RR: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end


%Storage installed less than capacity
StorageNodeK=ones(1,NumberBuses);
%LI-ION
b_vector=[b_vector;LIION_SC];
sense_vector=[sense_vector;'<'];
A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
model.constrnames{MatrixIndex} = sprintf('Installed LI-ION<Capacity');
MatrixIndex=MatrixIndex+1;

StorageNodeK=zeros(1,NumberBuses);

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
                A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_A_delta StorageNodeK];
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
                A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_A_delta StorageNodeK];
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
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;GeneratorCapacity(i)];
    end
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
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;LIION_DR];
    end
end
%LI-ION Discharge Rate
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;LIION_DR];
    end
end
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

model.lb=LowerBound;
model.ub=UpperBound;


%% Gurobi Optimization
tic
result=gurobi(model,params);
toc

%% Plot Results
ProductionBus1=result.x(1:TotalSimulation);
ProductionBus2=result.x(TotalSimulation+1:2*TotalSimulation);
ProductionBus3=result.x(2*TotalSimulation+1:3*TotalSimulation);
ProductionBus4=result.x(3*TotalSimulation+1:4*TotalSimulation);
ProductionBus5=result.x(4*TotalSimulation+1:5*TotalSimulation);
ProductionBus6=result.x(5*TotalSimulation+1:6*TotalSimulation);
ProductionBus7=result.x(6*TotalSimulation+1:7*TotalSimulation);
LiIONStorageLevelBus1=result.x(7*TotalSimulation+1:8*TotalSimulation);
LiIONStorageLevelBus2=result.x(8*TotalSimulation+1:9*TotalSimulation);
LiIONStorageLevelBus3=result.x(9*TotalSimulation+1:10*TotalSimulation);
LiIONStorageLevelBus4=result.x(10*TotalSimulation+1:11*TotalSimulation);
LiIONStorageLevelBus5=result.x(11*TotalSimulation+1:12*TotalSimulation);
LiIONStorageLevelBus6=result.x(12*TotalSimulation+1:13*TotalSimulation);
LiIONStorageLevelBus7=result.x(13*TotalSimulation+1:14*TotalSimulation);
LiIONChargeRateBus1 = result.x(14*TotalSimulation+1:15*TotalSimulation);
LiIONChargeRateBus2 = result.x(15*TotalSimulation+1:16*TotalSimulation);
LiIONChargeRateBus3 = result.x(16*TotalSimulation+1:17*TotalSimulation);
LiIONChargeRateBus4 = result.x(17*TotalSimulation+1:18*TotalSimulation);
LiIONChargeRateBus5 = result.x(18*TotalSimulation+1:19*TotalSimulation);
LiIONChargeRateBus6 = result.x(19*TotalSimulation+1:20*TotalSimulation);
LiIONChargeRateBus7 = result.x(20*TotalSimulation+1:21*TotalSimulation);
LiIONDisChargeRateBus1 = result.x(21*TotalSimulation+1:22*TotalSimulation);
LiIONDisChargeRateBus2 = result.x(22*TotalSimulation+1:23*TotalSimulation);
LiIONDisChargeRateBus3 = result.x(23*TotalSimulation+1:24*TotalSimulation);
LiIONDisChargeRateBus4 = result.x(24*TotalSimulation+1:25*TotalSimulation);
LiIONDisChargeRateBus5 = result.x(25*TotalSimulation+1:26*TotalSimulation);
LiIONDisChargeRateBus6 = result.x(26*TotalSimulation+1:27*TotalSimulation);
LiIONDisChargeRateBus7 = result.x(27*TotalSimulation+1:28*TotalSimulation);

DemandAvg=DemandAvg(1:TotalSimulation);
ExpArgAvg=ExpArgAvg(1:TotalSimulation);
ExpBrasilAvg=ExpBrasilAvg(1:TotalSimulation);

WindAvg=WindAvg(1:TotalSimulation);
SolarAvg=SolarAvg(1:TotalSimulation);
BiomassAvg=BiomassAvg(1:TotalSimulation);
ThermalAvg=ThermalAvg(1:TotalSimulation);

TotalProduction=ProductionBus1+ProductionBus2+ProductionBus3+ProductionBus4+ProductionBus5+ProductionBus6+ProductionBus7;
TotalStorageLevel=LiIONStorageLevelBus1+LiIONStorageLevelBus2+LiIONStorageLevelBus3+LiIONStorageLevelBus4+LiIONStorageLevelBus5+LiIONStorageLevelBus6+LiIONStorageLevelBus7;
TotalChargeRate=LiIONChargeRateBus1+LiIONChargeRateBus2+LiIONChargeRateBus3+LiIONChargeRateBus4+LiIONChargeRateBus5+LiIONChargeRateBus6+LiIONChargeRateBus7;
TotalDisChargeRate=LiIONDisChargeRateBus1+LiIONDisChargeRateBus2+LiIONDisChargeRateBus3+LiIONDisChargeRateBus4+LiIONDisChargeRateBus5+LiIONDisChargeRateBus6+LiIONDisChargeRateBus7;

%Plot Overall System
figure
plot(DemandAvg+ExpArgAvg+ExpBrasilAvg,'k--')
hold on
plot(WindAvg+SolarAvg+BiomassAvg,'k')
hold on
plot(TotalProduction)
hold on
plot(TotalStorageLevel)
hold on
plot(TotalChargeRate,'--')
hold on
plot(TotalDisChargeRate)
hold on
legend('Demand+Exportation','Wind+Solar+Biomass','Production','LiIONStorage','LiIONCharge','LiIONDisCharge')

%Plot Storage
figure
plot(TotalStorageLevel)
hold on
plot(TotalChargeRate,'--')
hold on
plot(TotalDisChargeRate)
hold on
legend('LiIONStorage','LiIONCharge','LiIONDisCharge')

%Plot Thermal
figure
plot(ThermalAvg)
hold on
plot(ProductionBus4+ProductionBus5)
legend('Thermal pre-battery','Thermal post-battery')

%Plot Error
figure
TotalDemand=DemandAvg+ExpArgAvg+ExpBrasilAvg+TotalChargeRate;
TotalPower=WindAvg+SolarAvg+BiomassAvg+TotalProduction+TotalDisChargeRate;
plot(TotalDemand-TotalPower)

fprintf('Li-ION Siting Bus1 = %.2f \n',result.x(35*TotalSimulation+1))
fprintf('Li-ION Siting Bus2 = %.2f \n',result.x(35*TotalSimulation+2))
fprintf('Li-ION Siting Bus3 = %.2f \n',result.x(35*TotalSimulation+3))
fprintf('Li-ION Siting Bus4 = %.2f \n',result.x(35*TotalSimulation+4))
fprintf('Li-ION Siting Bus5 = %.2f \n',result.x(35*TotalSimulation+5))
fprintf('Li-ION Siting Bus6 = %.2f \n',result.x(35*TotalSimulation+6))
fprintf('Li-ION Siting Bus7 = %.2f \n',result.x(35*TotalSimulation+7))
