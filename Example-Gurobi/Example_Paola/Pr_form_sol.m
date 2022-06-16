% 
% 
%  Defining and solving the Optimization problem for line reinforcement
%

% existing generators
gen_exi.ng=ng;
gen_exi.cg=cg;  
gen_exi.bgen=bgen;
gen_exi.Fuelc=Fuelc;
gen_exi.nrg=nrg;                     % Number of random generators
gen_exi.cgr=cgr;                     % Cost of the random generators
gen_exi.bgenr=bgenr;                 % Bus of random generators  




% scenario parameters
  em=stagenode(NE,Nem);   %each entry i of em contains the stage to which node i belongs to


  
% Network data 

net_data.nl=nl;
net_data.nb=nb;
net_data.xr=xr;
net_data.Ob=Ob;
net_data.Db=Db;
net_data.wind=wind;
net_data.dtn=dtn;

% line investment parameters

Fmax=kron(Fmxw,Line_candidate);          % Matrix Nw \times NL; Maximum capacity provided for each line and project F_{l,w}^max
glw=kron(FCT,Line_candidate.*Ll);        % Matrix Nw \times NL where each entry is the annuitized fixed investment 
                                          % cost for line l and option w    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Existing storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sto_exis.nse=nse;
sto_exis.bsto=bsto;
sto_exis.char_eff=char_effe;
sto_exis.disch_eff=disch_effe;
sto_exis.hbar_che=hbar_che;
sto_exis.hbar_dise=hbar_dise;
sto_exis.nubare=nubare;
sto_exis.disch_coste=disch_coste;
sto_exis.charg_coste=charg_coste;
sto_exis.flaghNTE=flaghNTE; 
sto_exis.flagh0E=flagh0E; 
sto_exis.h0E=h0E;


% computing the the sets of all parents nodes of m  for each stage and its
% probability

[p0m,pim]=phim(NM,em,fathernode,Nodeprob); 
   

% Compute the discount factors

YNE=NYE*NE;         % Total number of years considered in the project.

% The computation of the cumulative discount assumes that the number of
% years in each stage are identical.

ir=[0:YNE-1];
ye=[1:NYE:YNE];

ai=1./((1+r).^ir);  % each entry contain the discount for the year i;



for i=1:NE
   reI(1,i)=sum(ai((i-1)*NYE+1:YNE));                % Cumulative discount factor for line investment  cost in epoch e. The lenght of the vector is NE  
   reO(1,i)=sum(ai((i-1)*NYE+1:i*NYE));             % Cumulative discount factor for operation   cost in epoch e 
end


% Cumulative discount relative to each scenario node

reOm=reO(1,em(1:end));
reIm=reI(1,em(1:end));


scenario_data.p0m=p0m;
scenario_data.reOm=reOm;
scenario_data.reIm=reIm;
scenario_data.pim=pim;
scenario_data.NM=NM;
scenario_data.reO=reO;
scenario_data.Nem=Nem;
scenario_data.em=em;


[Alconst,object,bterm,ssense,vtype,lb,ub,ex_st_map]=define_probl(net_data,gen_exi,NMoc,taut,GG,NT,NBl,Fmax,Line_candidate,NWl,pmaxmn,glw,F0,NTb,sto_exis,scenario_data);
                                                                


model.A =Alconst;
model.obj = full(object);
model.modelsense = 'Min';
model.rhs = full(bterm);
model.sense=ssense;
model.vtype=vtype;
model.lb=full(lb);
model.ub=full(ub);




 start=tic;
 result = gurobi(model,params);
 extimeUn=toc(start);
 

 

 
 disp(result.objval);



 nvinv=nl*NWl*NM;
 

 


 
% Retrieve binary variables signifying expansion option w for line l in
% decision point m
   
  extracbeta=sparse(eye(nl));
  exbetam=kron(eye(NWl*NM),extracbeta);
  betamwl=exbetam*result.x(1:nvinv);
  betal=sparse(reshape(betamwl,nl,NWl*NM));
  full(betal)

% the variable betal contains information on line reinforcements (investment).
% It is a matrix  of dimension (nb \times NWl NM)
% If the number of investment options is two and NM is 3 the first two
% columns contains the possible investment at each bus node for the first
% scenario node m=1 and the two investment options. The columns 3 and 4 contains the investment at
% scenario node m=1 for the two investment options NWl and so on.




  
  
  
