
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Three bus data   (structure presented in Ioannis' thesis
%  (pages 122-123))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear model;
clear all


% Horizon parameters
 
NE=3;   % Number of stages Stages
NYE=5;  % Numbers of years in a stage. 


%  Setting Gurobi parameters

% params.outputflag=0;
% params.method=4; 
% params.intfeastol=1e-9;
% params.optimalitytol=1e-9;
params.TimeLimit=24*3600;  % Set execution time limit



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data on the installed generators
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ng=3;                             % number of installed generators
                                   % number of random (wind) generators
cg=[200 200 1500]./100;           % generator's capacity (MW)  (Since the line reactance is given in p.u. the MW have been divided by 100  ...)
bgen=[1 2 3];                     % generator's bus
Fuelc=[30 35 40].*100;            % Operation Cost £/MWh


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data on the random generators, related variables and scenarios
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nrg=1;                 % Bus of wind generator
cgr=0.05.*100;         % Wind generator's Fuel Cost £/MWh
NM=7;                  % number of nodes,  numbered sequentially 
bgenr=[3];             % generator's bus for the additional generator.
pmaxmn=[0 400 200 600 400 400 200]./100;   % evolution of the random varaiable (MW)
Nem=[1 2 4];                    % The i-entry of the vector Nem contains the number of scenarios tree nodes at epoch e; The lenght of Nem is NE;   
fathernode=[0,1,1,2,2,3,3];      %the father's node of node i is fathernode(i)

% Probability of each arch

Nodeprob=sparse(zeros(NM,NM));
Nodeprob(1,2)=0.6;   % Denote a connection from node 1 to node 2 with probability 0.6; 
Nodeprob(1,3)=0.4;
Nodeprob(2,4)=0.3;
Nodeprob(2,5)=0.7;
Nodeprob(3,6)=0.2;
Nodeprob(3,7)=0.8;

checkpb=sum(Nodeprob'); % add check error control

em=stagenode(NE,Nem);   %each entry i of em contains the stage to which node i belongs to

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Network Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nl=3;  % number of lines
nb=3;  % number of busses 

xr=[0.2 0.2 0.2]; % reactance of the line [p.u.]
Ll=[100 100 100]; % Line length (km) 
%Ic=[100 100 100]; % Initial capacity (MW)
Ob=[1 1 2];       % Origin bus (u_l)
Db=[2 3 3];       % Destination bus (v_l)
F0=[100, 100, 100]./100;    % Initial capacity (MW)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investment data on Line Reinforcement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data on candidate projects

Line_candidate=ones(1,nl);                   % Line candidates for an investement with entries 1. The length of the vector must be nl 

NWl=2;                                       % Number of candidate projetcs for line l
Fmxw=[150; 300]./100;                    %Maximum capacity provided for each line and project  F_{w}^max                     
FCT=[40000; 60000];                      % Annuitized  fixed cost (£/(Km yr))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nse=2;                             %  number of installed storage units
bsto=[2 3];                        %  storage's bus
char_effe=[0.9 ,0.8];      % Vector of dimension nse containing the charging efficency coefficients for each installed storage. 
disch_effe=[0.7, 0.75]; 
disch_coste=[5, 4 ].*100;       % Vector of dimension nse containing the discharging cost for each operative storage. (£/MW)
charg_coste=[3, 3.5].*100  ;  % Vector of dimension nse containing the charging cost for each operative storage. (£/MW)
hbar_che=[400, 500]./100;           % Maximum charge rate of storage device (MW)  for each storage techn
hbar_dise=[600, 450]./100;          % Maximum discharge rate of storage device (MW)  for each operative storage.
nubare=[1600,  2000]./100;       % Energy capacity of storage device (MWh) for each operative storage. 
h0E=(3/4)*nubare;
flaghNTE=[1 1];    % The entries of this vector can be set 1 or 0. If the entry of flaghNTE is one the initial condition h0 for the storage is assigned to be equal the terminal state hNT 
                  % If the entry of flaghNT=0 the terminal state is free and it has to satisfies the phisical limitations only. 
flagh0E=[0  0];    %this variable can be set 1 or 0. If flagh0E is one the initial condition for the storage is assigned to be equal h0 
                  % otherwise is a variable of the optimization problem for each storage technology

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Information on demand blocks, demand and wind data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


NBl=5; %Number of blocks
               

NTb=[24, 24, 24 24, 24];  % Number of samples in each demand block
NT=sum(NTb);    
staut=[1,1 1 1 1];         %sampling time (h) for each block period.

taut{1}= ones(1,NTb(1))*staut(1);       
taut{2}= ones(1,NTb(2))*staut(2);  
taut{3}= ones(1,NTb(3))*staut(3);  
taut{4}= ones(1,NTb(4))*staut(4);  
taut{5}= ones(1,NTb(4))*staut(4);

NMoc=[12.5 12.5 12.5 12.5 2];  % Number of occurrence of each block 


% load data on demand and wind 
  
load datademand                  

dtn=round(peakdemand.*kron(loadfactor,bus_percentage))./100;  % Demand in MW converted p.u.
% Matrix NT \times NB where each entry is the demand (p.u.) in period t at node (Bus) n


% wind coefficients a vector NT \times 1

% Note datademands already contains the vector with the information on wind


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Additional data 

GG=30000.*100;             % Value of lost load £/MWh 
r=5/100;                   % Interest rate (Usually 5% or 3%)




