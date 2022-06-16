%%%%%%%%%%%%%%%%%%
%% THREE BUS TEST%
%%%%%%%%%%%%%%%%%%

clear all
clc

%%%%%%%%%%%%%%%%
% NETWORK DATA %
%%%%%%%%%%%%%%%%

NumberBuses=3;  % number of buses 

xr=[0.05 0.05 0.05]; % reactance of the line [p.u.]
Ob=[1 1 2];       % Origin bus
Db=[2 3 3];       % Destination bus
F0=[1, 1, 1]./100;    % Lines capacity (MW)

%%%%%%%%%%%%%%%%%%
% GENERATOR DATA %
%%%%%%%%%%%%%%%%%%

NumberGenerators=3;                           % number of installed generators
GeneratorCapacity=[2000 2000 2000]./100;        % generator's capacity (MW)
BusGenerator=[1 2 3];                         % generator's bus
FuelCost=[100 100 100].*100;                     % Operation Cost £/MWh


%%%%%%%%%%%%%%%%%%%%%%%
% WIND GENERATOR DATA %
%%%%%%%%%%%%%%%%%%%%%%%

%Wind data over a day
load WindAndDemandData
WindAvg=WindAvg/100; %p.u

%Wind in each bus
WindBus=[0 0 1];   

%%%%%%%%%%%%%%%%
% STORAGE DATA %
%%%%%%%%%%%%%%%%

%Pump Storage Hydro
PSH_ef= 0.87;         % Efficiency
PSH_SC= 1500/100;         % Storage Capacity (MWh)
PSH_DR= 250/100;          % Dis/Charge Rate (MW)

%LI-ION Batteries
LIION_ef=0.94;
LIION_SC=500/100;
LIION_DR=50/100;

%%%%%%%%%%%%%%
% OTHER DATA %
%%%%%%%%%%%%%%

TimeStep=10/60;                        %hours
TotalSimulation=24/TimeStep;   %simlulate over a day 
RampRate=4/100;                        %MW/time step

%Demand data
DemandBus=[0.20 0.20 0.60];
DemandAvg=DemandAvg/100; %p.u

%Gurobi Parameters
params.TimeLimit=24*3600;  % Set execution time limit
% Build model
model.modelname  = 'Storage Siting';
model.modelsense = 'min';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FORMULATE THE PROBLEM %
%%%%%%%%%%%%%%%%%%%%%%%%%%

HourlyCost=[];
PumpHydroStorageVar=[];
LIIONStorageVar=[];
RateChargePHSVar=[];
RateDisChargePHSVar=[]; 
RateChargeLIIONVar=[]; 
RateDisChargeLIIONVar=[]; 
AngleDeltaVar=[];
StorageNodeK=zeros(1,2*NumberBuses);


%% Define Variables
%Production
index_var=1;
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('Production: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%PH Storage Level
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('PH Storage: Bus%d Time%d',k,i);
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

%PH Charge Rate
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('PH Charge Rate: Bus%d Time%d',k,i);
        index_var=index_var+1;
    end
end

%PH Discharge Rate
for k=1:NumberBuses
    for i=1:TotalSimulation
        model.varnames{index_var} = sprintf('PH Discharge Rate: Bus%d Time%d',k,i);
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

%K-siting PH
for k=1:NumberBuses
        model.varnames{index_var} = sprintf('K siting PH: Bus%d',k);
        index_var=index_var+1;
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
    
    %Pump Hydro
    PumpHydroStorage{k}=zeros(1,TotalSimulation);
    PumpHydroStorageVar=[PumpHydroStorageVar 0*PumpHydroStorage{k}];                %Variable for optimization problem
    
    %LI-ION Batteries
    LIIONStorage{k}=zeros(1,TotalSimulation);
    LIIONStorageVar=[LIIONStorageVar 0*PumpHydroStorage{k}];                        %Variable for optimization problem
    
    %Rate Charge PHS
    RateChargePHS{k}=zeros(1,TotalSimulation);
    RateChargePHSVar=[RateChargePHSVar 0*RateChargePHS{k}];                  %Variable for optimization problem
 
    %Rate DisCharge PHS
    RateDisChargePHS{k}=ones(1,TotalSimulation); 
    RateDisChargePHSVar=[RateDisChargePHSVar 0*RateDisChargePHS{k}];         %Variable for optimization problem
    
    %Rate Charge LIION
    RateChargeLIION{k}=zeros(1,TotalSimulation); 
    RateChargeLIIONVar=[RateChargeLIIONVar 0*RateChargeLIION{k}];           %Variable for optimization problem 
    
    %Rate DisCharge LIION
    RateDisChargeLIION{k}=zeros(1,TotalSimulation); 
    RateDisChargeLIIONVar=[RateDisChargeLIIONVar 0*RateDisChargeLIION{k}];  %Variable for optimization problem

    %Rate DisCharge LIION
    AngleDelta{k}=zeros(1,TotalSimulation); 
    AngleDeltaVar=[AngleDeltaVar 0*AngleDelta{k}];  %Variable for optimization problem
end

%Define Objective Function: Minimize Production Cost
%Rest of the variables for optimization constrains
model.obj=[HourlyCost PumpHydroStorageVar LIIONStorageVar RateChargePHSVar RateDisChargePHSVar RateChargeLIIONVar RateDisChargeLIIONVar AngleDeltaVar StorageNodeK];


%% Build A Matrix - Equality Equations

A=zeros((TotalSimulation-1)*2*NumberBuses,NumberBuses*TotalSimulation*8+length(StorageNodeK)); %8 comes from model.obj 

aux_zeros=zeros(1,NumberBuses*TotalSimulation);

%Initialize b vector (Ax=b) and sense vector (=,<,>)
b_vector=[];
sense_vector=[];

MatrixIndex=1;

%Storage level Equation: s(t)=s(t-1)+(ef_ch.rate_ch(t)-rate_dis(t)/ef_dis)deltaT
for n=1:NumberBuses
for k=1:(TotalSimulation-1)
    %Create constrains for Pump Hydro Storage
    aux_A_storage_PH=zeros(1,NumberBuses*TotalSimulation);
    aux_A_charge_PH=zeros(1,NumberBuses*TotalSimulation);
    aux_A_discharge_PH=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_PH(TotalSimulation*(n-1)+k)=1;
    aux_A_storage_PH(TotalSimulation*(n-1)+k+1)=-1;
    aux_A_charge_PH(TotalSimulation*(n-1)+k+1)=PSH_ef*TimeStep;
    aux_A_discharge_PH(TotalSimulation*(n-1)+k+1)=-TimeStep/PSH_ef;
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_PH aux_zeros aux_A_charge_PH aux_A_discharge_PH aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('PH Storage Level: Bus%d Time%d', n, k+1);
    %Create constrains for batteries storage
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_charge_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_discharge_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    aux_A_storage_LIION(TotalSimulation*(n-1)+k+1)=-1;
    aux_A_charge_LIION(TotalSimulation*(n-1)+k+1)=LIION_ef*TimeStep;
    aux_A_discharge_LIION(TotalSimulation*(n-1)+k+1)=-TimeStep/LIION_ef;
    A(MatrixIndex+(TotalSimulation-1)*NumberBuses,:)=[aux_zeros aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_A_charge_LIION aux_A_discharge_LIION aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex+(TotalSimulation-1)*NumberBuses} = sprintf('LI-ION Storage Level: Bus%d Time%d', n, k+1);
    b_vector=[b_vector;0];
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end
end
MatrixIndex=MatrixIndex+(TotalSimulation-1)*NumberBuses; %Because of two storage technologies
%Cyclic storage level Equation: s(t=1)=s(t=T)
%Pump Hydro
for n=1:NumberBuses
    aux_A_storage_PH=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_PH(TotalSimulation*(n-1)+1)=1;
    aux_A_storage_PH(TotalSimulation*(n-1)+TotalSimulation)=-1;
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_PH aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Cyclic PH Storage: Bus%d', n);
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
    MatrixIndex=MatrixIndex+1;
end
%Battery
for n=1:NumberBuses
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+1)=1;
    aux_A_storage_LIION(TotalSimulation*(n-1)+TotalSimulation)=-1;
    A(MatrixIndex,:)=[aux_zeros aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
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
            
        %RATE CHARGE-DISCHARGE PHS
        aux_A_charge_PH=zeros(1,NumberBuses*TotalSimulation);
        aux_A_discharge_PH=zeros(1,NumberBuses*TotalSimulation);
        aux_A_charge_PH(TotalSimulation*(i-1)+m)=-1;
        aux_A_discharge_PH(TotalSimulation*(i-1)+m)=1;
        
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
        
        %WIND(t)
        Wind=WindAvg(m)*WindBus(i);
        end
        
    A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_A_charge_PH aux_A_discharge_PH aux_A_charge_LIION aux_A_discharge_LIION aux_A_delta StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Power Balance: Bus%d Time%d', i, m);
    MatrixIndex=MatrixIndex+1;  
    b_vector=[b_vector;Demand-Wind];
    sense_vector=[sense_vector;'='];
    end
end

%Finally set slack bus angle = 0
for i=1:TotalSimulation
    [maxCapacity, SwingBus]= max(GeneratorCapacity);
    aux_A_delta=zeros(1,NumberBuses*TotalSimulation);
    aux_A_delta(TotalSimulation*(SwingBus-1)+i)=1;
    A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_A_delta StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Swing Angle: Time%d', i);
    MatrixIndex=MatrixIndex+1;  
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'='];
end


%% A Matrix - Inequality Equations

%Storage level (t) > 0 Pump Hydro
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for Pump Hydro Storage
    aux_A_storage_PH=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_PH(TotalSimulation*(n-1)+k)=1;
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_PH aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Storage Level PH > 0: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end

%Storage level (t) > 0 LI-ION
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for LI-ION Storage
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[aux_zeros aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Storage Level LI-ION >0: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end

%Storage level (t) < Kj Pump Hydro
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for Pump Hydro
    aux_A_storage_PH=zeros(1,NumberBuses*TotalSimulation);
    StorageNodeK=zeros(1,2*NumberBuses);
    aux_A_storage_PH(TotalSimulation*(n-1)+k)=1;
    StorageNodeK(n)=-1;
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[aux_zeros aux_A_storage_PH aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Storage Level PH<kj: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end

%Storage level (t) < Kj LI-ION
for n=1:NumberBuses
for k=1:TotalSimulation
    %Create constrains for LI-ION
    aux_A_storage_LIION=zeros(1,NumberBuses*TotalSimulation);
    StorageNodeK=zeros(1,2*NumberBuses);
    aux_A_storage_LIION(TotalSimulation*(n-1)+k)=1;
    StorageNodeK(NumberBuses+n)=-1; 
    b_vector=[b_vector;0];
    sense_vector=[sense_vector;'<'];
    A(MatrixIndex,:)=[aux_zeros aux_zeros aux_A_storage_LIION aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('Storage Level LI-ION <kj: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end
StorageNodeK=zeros(1,2*NumberBuses);

%Ramp Rate Generators > -RR
for n=1:NumberBuses
for k=1:TotalSimulation-1
    aux_A_production=zeros(1,NumberBuses*TotalSimulation);
    aux_A_production(TotalSimulation*(n-1)+k)=-1;
    aux_A_production(TotalSimulation*(n-1)+k+1)=1;
    b_vector=[b_vector;-RampRate];
    sense_vector=[sense_vector;'>'];
    A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
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
    A(MatrixIndex,:)=[aux_A_production aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
    model.constrnames{MatrixIndex} = sprintf('p(t)-p(t-1)<RR: Bus%d Time%d', n, k);
    MatrixIndex=MatrixIndex+1;
end
end


%Storage installed less than capacity
StorageNodeK=zeros(1,2*NumberBuses);
%Pump Hydro
for n=1:NumberBuses
    StorageNodeK(n)=1;
end
b_vector=[b_vector;PSH_SC];
sense_vector=[sense_vector;'<'];
A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
model.constrnames{MatrixIndex} = sprintf('Installed PH<Capacity');
MatrixIndex=MatrixIndex+1;

StorageNodeK=zeros(1,2*NumberBuses);

%LI-ION
for n=1:NumberBuses
    StorageNodeK(NumberBuses+n)=1;
end
b_vector=[b_vector;LIION_SC];
sense_vector=[sense_vector;'<'];
A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros StorageNodeK];
model.constrnames{MatrixIndex} = sprintf('Installed LI-ION<Capacity');
MatrixIndex=MatrixIndex+1;

StorageNodeK=zeros(1,2*NumberBuses);

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
                A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_A_delta StorageNodeK];
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
                A(MatrixIndex,:)=[aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_zeros aux_A_delta StorageNodeK];
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
%PH Storage Level
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;inf];
    end
end
%LI-ION Storage Leve
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;inf];
    end
end
%PH Charge Rate
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;PSH_DR];
    end
end
%PH Discharge Rate
for i=1:NumberBuses
    for k=1:TotalSimulation
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;PSH_DR];
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
%K-siting PH
for i=1:NumberBuses
    LowerBound=[LowerBound;0];
    UpperBound=[UpperBound;inf];
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
PHStorageLevelBus1=result.x(3*TotalSimulation+1:4*TotalSimulation);
PHStorageLevelBus2=result.x(4*TotalSimulation+1:5*TotalSimulation);
PHStorageLevelBus3=result.x(5*TotalSimulation+1:6*TotalSimulation);
LiIONStorageLevelBus1=result.x(6*TotalSimulation+1:7*TotalSimulation);
LiIONStorageLevelBus2=result.x(7*TotalSimulation+1:8*TotalSimulation);
LiIONStorageLevelBus3=result.x(8*TotalSimulation+1:9*TotalSimulation);
PHChargeRateBus1 = result.x(9*TotalSimulation+1:10*TotalSimulation);
PHChargeRateBus2 = result.x(10*TotalSimulation+1:11*TotalSimulation);
PHChargeRateBus3 = result.x(11*TotalSimulation+1:12*TotalSimulation);
PHDisChargeRateBus1 = result.x(12*TotalSimulation+1:13*TotalSimulation);
PHDisChargeRateBus2 = result.x(13*TotalSimulation+1:14*TotalSimulation);
PHDisChargeRateBus3 = result.x(14*TotalSimulation+1:15*TotalSimulation);
LiIONChargeRateBus1 = result.x(15*TotalSimulation+1:16*TotalSimulation);
LiIONChargeRateBus2 = result.x(16*TotalSimulation+1:17*TotalSimulation);
LiIONChargeRateBus3 = result.x(17*TotalSimulation+1:18*TotalSimulation);
LiIONDisChargeRateBus1 = result.x(18*TotalSimulation+1:19*TotalSimulation);
LiIONDisChargeRateBus2 = result.x(19*TotalSimulation+1:20*TotalSimulation);
LiIONDisChargeRateBus3 = result.x(20*TotalSimulation+1:21*TotalSimulation);


%For Testing
DemandAvg=DemandAvg(1:TotalSimulation);
WindAvg=WindAvg(1:TotalSimulation);

%Plots
figure
plot(DemandAvg,'k--')
hold on
plot(WindAvg,'k')
hold on
plot(ProductionBus1+ProductionBus2+ProductionBus3)
hold on
plot(PHStorageLevelBus1+PHStorageLevelBus2+PHStorageLevelBus3)
hold on
plot(LiIONStorageLevelBus1+LiIONStorageLevelBus2+LiIONStorageLevelBus3)
hold on
plot(PHChargeRateBus1+PHChargeRateBus2+PHChargeRateBus3,'--')
hold on
plot(PHDisChargeRateBus1+PHDisChargeRateBus2+PHDisChargeRateBus3)
hold on
plot(LiIONChargeRateBus1+LiIONChargeRateBus2+LiIONChargeRateBus3,'--')
hold on
plot(LiIONDisChargeRateBus1+LiIONDisChargeRateBus2+LiIONDisChargeRateBus3)
hold on
legend('Demand','Wind','Production','PHStorage','LiIONStorage','PHCharge','PHDischarge','LiIONCharge','LiIONDisCharge')


%Plots
figure
plot(PHStorageLevelBus1+PHStorageLevelBus2+PHStorageLevelBus3)
hold on
plot(PHChargeRateBus1+PHChargeRateBus2+PHChargeRateBus3,'--')
hold on
plot(PHDisChargeRateBus1+PHDisChargeRateBus2+PHDisChargeRateBus3)
legend('PHStorage','PHCharge','PHDischarge')

figure
plot(LiIONStorageLevelBus1+LiIONStorageLevelBus2+LiIONStorageLevelBus3)
hold on
plot(LiIONChargeRateBus1+LiIONChargeRateBus2+LiIONChargeRateBus3,'--')
hold on
plot(LiIONDisChargeRateBus1+LiIONDisChargeRateBus2+LiIONDisChargeRateBus3)
hold on
legend('LiIONStorage','LiIONCharge','LiIONDisCharge')

TotalDemand=DemandAvg'+PHChargeRateBus1+PHChargeRateBus2+PHChargeRateBus3+LiIONChargeRateBus1+LiIONChargeRateBus2+LiIONChargeRateBus3;
TotalPower=WindAvg'+ProductionBus1+ProductionBus2+ProductionBus3+PHDisChargeRateBus1+PHDisChargeRateBus2+PHDisChargeRateBus3+LiIONDisChargeRateBus1+LiIONDisChargeRateBus2+LiIONDisChargeRateBus3;
figure
plot(TotalDemand-TotalPower)

fprintf('PH Siting Bus1 = %.2f \n',result.x(24*TotalSimulation+1))
fprintf('PH Siting Bus2 = %.2f \n',result.x(24*TotalSimulation+2))
fprintf('PH Siting Bus3 = %.2f \n',result.x(24*TotalSimulation+3))
fprintf('Li-ION Siting Bus1 = %.2f \n',result.x(24*TotalSimulation+4))
fprintf('Li-ION Siting Bus2 = %.2f \n',result.x(24*TotalSimulation+5))
fprintf('Li-ION Siting Bus3 = %.2f \n',result.x(24*TotalSimulation+6))
