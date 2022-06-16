

function [Alconst,object,bterm,ssense,vtype,lb,ub,ex_st_map]=define_probl(net_data,gen_exi,NMoc,taut,GG,NT,NBl,Fmax,Line_candidate,NWl,pmaxmn,glw,F0,NTb,sto_exis,scenario_data)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute constraints for every scenario node m and  demand period t. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Predefintion of useful variables

Aeq='=';
Aleq='<';
vab='B';   % it is used to define a variable as Binary
vac='C';   % it is used to define a variable as continuous



%Assumption: a generator is directly connected to a single bus

% Vector of variables for a given t and a given node is
%   x_{m^*,t*}=[p_g(1),...,p_g(ngTot),f_l(1),......f_l(nl),\theta(2),....
% \theta(nb),d^*(1), .....,d^*(nb)]
% The angle at the first bus is set \theta(1)=0 and so it is not a variable
% of the optimization problem part of the cost related to the constraint
% cost and loss Load.
%
%
% So for each block we have a variables' vector defined by
%
%  OP_TBi[m^*] =[x_{m^*,t1},x_{m^*,t2}, ..., ,x_{m^*,TN}], 
%
%  OP_TM=[OP_TBi[1], .... ,    OP_TBi[NM]]                       
%
%  OP_Tot=[OP_TM,Stor_exist]
%   
% where   Stor_exist = [      ]  (see later)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALL THE OPERATION VARIABLES RELATED TO THE
% STORAGE ARE AT THE END OF THE VECTOR OF VARIABLES.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve information on existing generetors 

ng=gen_exi.ng;     
cg=gen_exi.cg;     
Fuelc=gen_exi.Fuelc;
bgen=gen_exi.bgen;
nrg=gen_exi.nrg;
cgr=gen_exi.cgr;
bgenr=gen_exi.bgenr;

      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  Retrieve indormation on existing storage 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nse=sto_exis.nse;
%cse=sto_exis.cse;
bsto=sto_exis.bsto;
char_effe=sto_exis.char_eff;
disch_effe=sto_exis.disch_eff;
hbar_che=sto_exis.hbar_che;
hbar_dise=sto_exis.hbar_dise;
nubare=sto_exis.nubare;
%helife=sto_exis.helife;
disch_coste=sto_exis.disch_coste;
charg_coste=sto_exis.charg_coste;
flaghNTE=sto_exis.flaghNTE; 
flagh0E=sto_exis.flagh0E; 
h0E=sto_exis.h0E;
%ahE=sto_exis.ahE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%  Retrieve indormation on Network data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nl=net_data.nl;
nb=net_data.nb;
xr=net_data.xr;
Ob=net_data.Ob;
Db=net_data.Db;
wind=net_data.wind;
dtn=net_data.dtn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Retrieve scenario data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reOm=scenario_data.reOm;
reIm=scenario_data.reIm;
pim=scenario_data.pim;
NM=scenario_data.NM;
reO=scenario_data.reO;
Nem=scenario_data.Nem;
em=scenario_data.em;
p0m=scenario_data.p0m;

% demand periods and nodes 

InM=incidmatrix(nb,nl,xr,Ob,Db);   %% computation of the incidence matrix to define the distribution of the
                                    % power in the network

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Defining useful matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
InoT=ones(1,NT);
InT=sparse(eye(NT));
Imi=sparse(eye(NM));

nbbr=nl+2*nb-1;                % Number of variables describing the network involved in power flow balance
% Power distribution term over the network according to the DCOPF formulation 
% There are nl constraints and nb-1+nl variables 

nrgtot=ng+nrg;            % number of existing  conventional generators and random wind generators

bgentot=[bgen,bgenr];     % bus number of existing conventional genrators and random wind generators
      
Fuelct=[Fuelc,cgr];      % Marginal cost 


noste=NM*(3*NT*nse+nse*NBl); % Total number of variables modelling existing storage  



%nocb=NT*(NM*(nbbr+nrg+nbpv+nWgn_Iw+nQBc)+3*Nem*gngE'+Nem*ngesj');    


Ntcong= nrgtot+nbbr;  % total number of variables describing line trasmission and generators at operational level for a given time t and scenario node m 

% creating a matrix filtering the operative (they must be operative and not retired) generators

  
% total number of  variables at operational level 

Noplg=NM*NT*Ntcong;
Nopvar=Noplg+noste;

% vtypeop=repmat(char(0),1,Noplg); % vector initialization
% objectng=sparse(1,Noplg);
objectbs=sparse(1,NT*Ntcong);


% indEst=0;
% jobj=0;
% 
% stepibe=0;
% 
% for ibe=1:NE 
%  %vibe=ones(1,ibe);
%  factNE=reO(ibe)*pim(indEst+1:indEst+Nem(ibe));
%  indEst= indEst+Nem(ibe);
% end 
 
stepng=0;
 for ibl=1:NBl
  % define part of cost at operational level    
  %objectbsa=NMoc(ibl)*kron(taut{ibl},[Fuelces(iges(indes{ibe})), Fuelc(ccgE{ibe}),cgr, cpv, Mwgcost,OgcE, sparse(1,nQBc), sparse(1,nl+nb-1),GG*ones(1,nb),zeros(1,nres(ibe)+nresEs(ibe)), Fuelc(ccgE{ibe}).*PgminE(ccgE{ibe})+Noloadc(ccgE{ibe}), zeros(1,2*gngE(ibe)), PmE{ibe}.*OgcE+NLcE,zeros(1,2*gennwm(ibe))]');

  objectbs(stepng+1:stepng+NTb(ibl)*Ntcong)=NMoc(ibl)*kron(taut{ibl},[Fuelct, sparse(1,nl+nb-1), GG*ones(1,nb)]);  
 
%   spptype=repmat([repmat(vac,1,Ntcong(ibe)),repmat(vab,1,Ntbint(ibe)+Ntbinex(ibe))],1,NTb(ibl));
%   vtypeo(1,stepng+1:stepng+NTb(ibl)*levng)=spptype;

  stepng=stepng+NTb(ibl)*Ntcong;

 end


% part of the cost related to the constraint cost
%constcost=sparse([kron(reOm.*pim, constcost), kron(reOm.*pim,zeros(1,nbstot*(2*NT+NBl)))]);


% Define the variable type for x_{m^*,t*}



% % Power flow over the line l which has to be equated to zero. 
 
Pfb=sparse([zeros(nl,nrgtot), -diag(ones(1,nl)), InM, zeros(nl,nb)]);
% 
% % define the constraints of the power flow over the line l for all demand
% % period t and all node m
% 

Pfbt=sparse(kron(eye(NT*NM),Pfb));


ssense=repmat(Aeq,nl*NT*NM,1);
bterm=sparse(nl*NT*NM,1);



% Lower and upper bounds on the variables modelling investments on line
% reinforcements

pbo=reshape(repmat(Line_candidate,NWl,1)',nl*NWl,1);

% Lower and upper bounds on the variables modelling operation costs



StepNTi=0;
upbt=sparse(zeros(1,Noplg));
lbu=sparse(zeros(Noplg,1));

for ij=1:NM 
      
     
    % upper bound
        
     upb=[cg, pmaxmn(:,ij)', inf*ones(1,nl), pi*ones(1,nb-1),inf*ones(1,nb)]; % upper bound for a given t for conventional generators and network variables
     upbr=repmat(upb,NT,1);
     upbr(:,ng+1:nrgtot)=upbr(:,ng+1:nrgtot).*[kron(ones(1,nrg),wind)];   
     
   
     upbt(1,StepNTi+1:StepNTi+NT*Ntcong)=reshape(upbr',1,NT*Ntcong);  % upper bound for all times 
    
     % lower bound for conventional generators and network variables
  
     lbu(StepNTi+1:StepNTi+NT*Ntcong,1)=[kron(InoT',[zeros(nrgtot,1);  -inf*ones(nl,1); -pi*ones(nb-1,1); 0*ones(nb,1)])];     
     StepNTi= StepNTi+NT*Ntcong;
   
end  

%ub=upbt';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  bounds on the storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% bounds on the existing storage
ubhe=sparse(noste,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Definition of cost and  constraints for the existing storage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Bus-to-storage incidence matrix of existing storage  


Nsm=nse*(3*NT+NBl);
oneg=ones(1,nse);
Sg=1:nse;
Sng=sparse(bsto,Sg,oneg,nb,nse,nse);


%Variables' initialization
% SE_{m^*,t^*}=[\tilde{h}_1,...,\tilde{h}_{nse},h^d_{1},...,h^d_{nse},h^c_{1},...,h^c_{nse}];
% shE_{m^*}=[SE_{m^*,tb_1},...,
% SE_{m^*,NTb(1)},SE_{m^*,tb_2},...,SE_{m^*,NTb(2)},....]

costeh=sparse(1,Nsm);
ubhem=sparse(Nsm,1);
bassel=[nubare, hbar_dise, hbar_che]';
baseceh=[zeros(1,nse), disch_coste, charg_coste];
baseKibS=[zeros(nb,nse), Sng, -Sng];
B_E=[diag(ones(1,nse)), zeros(nse, 2*nse)];         %  [\tilde{h}, h_d, h_c]; B contains a base for the state variables 
                                                    %   of the existing storage
ubasec=[sparse(nse,2*nse),diag(ones(1,nse))];                     % base for the charging inputs h_c
ubased=[sparse(nse,nse),diag(ones(1,nse)), sparse(nse,nse)];      % base for the discharging inputs h_d                                                   
                                                    

                                               
ubaseE=[zeros(nse,nse),-diag(ones(1,nse)./disch_effe), diag(char_effe.*ones(1,nse))];
chE=[];
KirbalSE=[];
h0hNT_E=[];

h0cE=[];
rnse=0;

% maps extracting the state of the storage  and the charging and discharging inputs   
ex_st_map.charge=sparse(NT*nse,3*NT*nse+nse*NBl);
ex_st_map.discharge=sparse(NT*nse,3*NT*nse+nse*NBl);
ex_st_map.sthtot=sparse(NT*nse+nse*NBl,3*NT*nse+nse*NBl);

hust=0;
hxst=0;

for ih=1:NBl

Ihb=eye(NTb(ih));
    
% Variables' upper bounds       
ubhem(rnse+1:rnse+nse*(3*NTb(ih)+1),1)=[nubare';repmat(bassel,NTb(ih),1)];

%Cost
costeh(1,rnse+nse+1:rnse+nse+nse*3*NTb(ih))=NMoc(ih)*[repmat(baseceh,1,NTb(ih))];


% Contribution of the existing storage to the balance equation

KirbalSE=blkdiag(KirbalSE,[repmat(zeros(nb,nse),NTb(ih),1), kron(eye(NTb(ih)),baseKibS)]);
 
% state equation constraints

B1E=kron(eye(NTb(ih)-1),B_E);
xp1=[zeros(NTb(ih)*nse,nse), kron(eye(NTb(ih)),B_E)];  
xp=[-blkdiag(diag(ones(1,nse)),B1E), zeros(NTb(ih)*nse,3*nse)];

uterm=[zeros(NTb(ih)*nse,nse),kron(diag(taut{ih}),ubaseE)];
 % ch is a matrix describing  the storage state equation
      
chE=blkdiag(chE,xp1+xp-uterm);

% boundary constraints  on the storage    
  h0hNT_E=blkdiag(h0hNT_E,[diag(ones(1,nse).*flaghNTE), sparse(nse,3*(NTb(ih)-1)*nse),-diag(ones(1,nse).*flaghNTE), sparse(nse,2*nse)]); 
  h0cE=blkdiag(h0cE,[diag(ones(1,nse).*flagh0E), sparse(nse,3*NTb(ih)*nse)]);
       
 
  % this part select the variables h^-_{m^*,t^*,n} and h^+_{m^*,t^*,n}
  % Selection matrix   
   ex_st_map.charge(hust+1:hust+nse*NTb(ih),rnse+1:rnse+nse*(3*NTb(ih)+1))=[sparse(nse*NTb(ih),nse), kron(Ihb,ubasec)];
   ex_st_map.discharge(hust+1:hust+nse*NTb(ih),rnse+1:rnse+nse*(3*NTb(ih)+1))=[zeros(nse*NTb(ih),nse), kron(Ihb,ubased)];
   ex_st_map.sthtot(hxst+1:hxst+nse*NTb(ih)+nse,rnse+1:rnse+nse*(3*NTb(ih)+1))=[eye(nse), zeros(nse,3*NTb(ih)*nse); zeros(nse*NTb(ih),nse), kron(Ihb,B_E)];
  
  
hust=hust+nse*NTb(ih);
hxst=hxst+nse*NTb(ih)+nse;  

% index update
rnse=rnse+nse*(3*NTb(ih)+1);

end


h0GE=repmat((h0E.*flagh0E)',NBl,1);  % right term of the costraints on the boundary conditions h0c of the existing storage


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% total number of investment variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvinv=nl*NWl*NM;


% total lower and upper bound

lb=[zeros(nvinv,1);lbu;sparse(noste,1)];%
ub=[repmat(pbo,NM,1); upbt';kron(ones(NM,1),ubhem)];

%*****************************



% % % Bus-to-generation incidence matrix of generators candidates
% % 

oneg=ones(1,nrgtot);
Jg=1:nrgtot;
Bng=sparse(bgentot,Jg,oneg,nb,nrgtot,nrgtot);

% Bus-to-line incidence matrix

% Inl=1  if the receiving bus for line l is n (row)
% Inl=-1 if the sending bus for line l is n

nlone=ones(1,nl);

btol=[-nlone, nlone];
nz_max=length(btol);
nlc=1:nl;
Ir=[nlc,nlc];
Ic=[Ob,Db];
Inl=(sparse(Ir,Ic,btol,nl,nb,nz_max))';

% Curtailed demand d^*_{m^*,t^*}(n) for each bus node 


% Kirbalan contains the system balance equation that must be equated to the
% demand at node n in period t^* (d_{t,n})

 Kirbalan=[Bng, Inl, zeros(nb,nb-1), diag(ones(1,nb))];
 Kirbalant=kron(eye(NT*NM),Kirbalan);


% Addition of storage contribution  to the balance equation 
 
 Kirbal=[Kirbalant,kron(Imi,KirbalSE)];
 
 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ssense=[ssense;repmat(Aeq,nb*NT*NM,1)];
dn=sparse(reshape(dtn',NT*nb,1));
bterm=[bterm;repmat(dn,NM,1)];

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
 
% Part defining complicating constraints in  "Complconstrt":  this part select the variable f_{m^*,t^*,l}

fmtlsc=sparse([diag(ones(1,nl))]);

%complicating constraints 

Complconstrb=sparse([zeros(nl,nrgtot),fmtlsc, zeros(nl,2*nb-1)]);
Complconstrt=kron(eye(NT*NM),Complconstrb);


sel=repmat(Aleq,nl*NT*NM,1);
seg=repmat(Aleq,nl*NT*NM,1);
termc=repmat(F0',NT*NM,1);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Investment modelling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% vectors of variables
% Line investment
%  
%  yInv_{m^*,w^*}=[\beta_1, ..., \beta_l];
%

%nvinv=nl*NWl*NM;   %line investment 

vtyl=repmat(vab,1,nvinv);
vtyht=repmat(vac,1,noste); % operational variables related to storage
vtypeop=repmat(vac,1,Noplg);

vtype=[vtyl,vtypeop, vtyht];


% first part of the cost
% Cost of the investment part

costinc=reshape(glw',1,NWl*nl);
winvs=reIm.*pim;

objecta=kron(winvs,costinc);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% constraint on the fixed cost related to the variables \beta_{m,l,w} for
% all m
%
%  \sum_{w} \beta(m^*)_{l,w}<=1
%


% Inizialization of the matrix describing the aggregate line capacity

cagg=nl*NWl;

% 
% exclusivity of line investment options

bmlw=diag(ones(1,nl));
bml=repmat(bmlw,1,NWl);
bmlM=kron(Imi,bml);

noteb=[ones(nl*NM,1)];
ssense=[ssense;repmat(Aleq,nl*NM,1)];


%%%%%%%%%%%%%%%%%
% Impose the structure for the aggragated capacity added up to node m
% to line l
%  Fmt contains the aggregate capacity \tilde\Beta^{inv}_{m,l} for an assigned m.  A different delay
%  for each line and project is considered.




%
% define the part of the complicated constaints relative to the investments
% \Beta^{inv}_{m,l}
Binvkl=sparse(repmat(bmlw,NT,1));

%%%%%%%%%%%%%
  % Finvmt is a matrix containing the information on   the aggregate
  % capacity F^{inv}_{m,l} for all nodes m and lines l 
 
Finvmt=sparse(NM*nl*NT,NM*nl*NWl);
findex=0;

for i=1:NM


pakmlw=p0m{i};
nnk=size(pakmlw,2);
 for l=1:nl
  for  jw=1:NWl
   for imk=1:nnk
     Finvmt(findex+1:findex+nl*NT,(pakmlw(imk)-1)*cagg+l+nl*(jw-1))=Fmax(jw,l)*Binvkl(:,l);
    end    
   end
  end
  
findex=findex+nl*NT;
end 


ssense=[ssense; sel;seg];
bterm=[bterm;noteb;termc;termc]; 


ssense=[ssense;repmat(Aeq,(nse*NT+2*nse*NBl)*NM,1)];
bterm=[bterm;zeros(nse*NT*NM,1);repmat([zeros(nse*NBl,1);h0GE],NM,1)];



% Total cost vector
objectb=sparse([kron(reOm.*pim,objectbs),kron(reOm.*pim,costeh)]);
object=[objecta,objectb];


%  define the matrix with the global constraints



Alconst=[sparse(NT*NM*nl,nvinv),   Pfbt, sparse(nl*NT*NM,noste); ...
         sparse(NT*NM*nb,nvinv),   Kirbal; ...
         bmlM, sparse(NM*nl,Nopvar); ...
         -Finvmt,  Complconstrt, sparse(nl*NT*NM,noste); ...
         -Finvmt, -Complconstrt, sparse(nl*NT*NM,noste); ...
          sparse(NT*nse*NM,nvinv+Noplg),kron(Imi,chE);...        % heqE
         sparse(2*nse*NM*NBl,nvinv+Noplg),  kron(Imi,[h0hNT_E;h0cE])];   %heqbE
          

     
      
     
      
   
      