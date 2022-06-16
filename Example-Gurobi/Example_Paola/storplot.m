%  evolution of candidate storage Deterministic 
%
%

% Select scenario node

sms=3;

% %select technology
% 
% tchs=1;
% 
% % select storage number nbs(tchs)
% 
% ssb=1;

% select time block
%
tbk=1

% start pointfor candidate technology



% parameters

nostem=(sms-1)*(3*NT*nse+nse*NBl);
nostemp=(3*NT*nse+nse*NBl);

% extract variables related to the selected options
% exinpst,exstor


gendl=NT*NM*(ng+nrg+nl+2*nb-1);
ins=nvinv+gendl+nostem;


varst=ex_st_map.sthtot*result.x(ins+1:ins+nostemp);      % extract the state of the storage devises   
varc=ex_st_map.charge*result.x(ins+1:ins+nostemp);       % extract the charging inputs 
varsd=ex_st_map.discharge*result.x(ins+1:ins+nostemp);   % extract the discharging inputs


stepNMu=NT*nse;
stepNMx=(NT+NBl)*nse;

Tfst=sum(NTb(1:tbk-1))*nse;

%indices to determine the beginning of the block tbk

jsindx=nse*(sum(NTb(1:tbk-1))+tbk-1);
jsindu=nse*sum(NTb(1:tbk-1));

for js=1:nse
    
stateh(:,js)=varst(jsindx+js:nse:jsindx+nse*(NTb(tbk)+1));
inputc(:,js)=varc(jsindu+js:nse:jsindu+nse*(NTb(tbk)));
inputd(:,js)=varsd(jsindu+js:nse:jsindu+nse*(NTb(tbk)));

end

timeu=[staut(tbk):staut(tbk):staut(tbk)*NTb(tbk)];
timex=[staut(tbk):staut(tbk):staut(tbk)*NTb(tbk)+staut(tbk)];




figure(1)
stairs((inputd-inputc).*100,'r','LineWidth',2)
xlabel('time')
ylabel('h_c-h_d,  MW')


figure(2)
plot(timex,stateh.*100)
hold on
stairs(timeu,(kron(char_effe,ones(NTb(tbk),1)).*inputc-inputd./(kron(disch_effe,ones(NTb(tbk),1)))).*100,'r')
%stairs((inputc-inputd).*100,'k--')
grid on
xlabel('time')
ylabel('MWh')

