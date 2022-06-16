% computation of the incidence matrix

function InM=incidmatrix(nb,nl,xr,Ob,Db)

% nl number of lines
% nb number of busses 
%xr reactance of the line [p.u.]     
% Ob Origin bus (u_l)
% Db Destination bus (v_l)

if nl>0

Tbn=nb-1;  % number of bus angles for n-1 nodes. The bus angle for the first node is set to 0


ObT=Ob-1;
DbT=Db-1;

nlO=find(ObT);
nlD=find(DbT);
lennl0=size(nlO,2);
lennlD=size(nlD,2);

slink=[ones(1,lennl0)./xr(nlO), -ones(1,lennlD)./xr(nlD)];
nz_max=length(slink);
Ir=[nlO,nlD];
Ic=[ObT(nlO),DbT(nlD)];



%S = sparse(I,J,S,m,n,maxnz). This creates an mxn sparse matrix with entry (I(k),J(k)) equal to Sk . 
%The optional argument maxnz causes Matlab to pre-allocate storage for maxnz nonzero entries, which can increase efficiency in the case when more nonzeros will be added later to S.

InM=sparse(Ir,Ic,slink,nl,Tbn,nz_max);

else
 InM=sparse(nl,nl);   
     
end
    