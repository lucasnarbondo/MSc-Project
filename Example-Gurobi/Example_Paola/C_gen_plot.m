% 
%   Plot - output of conventional generators and random plants 
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
tbk=1;

% start point for generators' ouputs

% parameters


% extract variables related to the selected options
Stepgd=(ng+nrg+nl+2*nb-1);
gendlM=NT*(sms-1)*Stepgd;
in_curt=(ng+nrg+nl+nb-1);

%gendl=NT*NM*(ng+nrg+nl+2*nb-1);

Tsumbk=sum(NTb(1:tbk-1));
ins=nvinv+gendlM+Stepgd*Tsumbk;


varopgt=result.x(ins+1:ins+NTb(tbk)*Stepgd);      

currt_dem=zeros(NTb(tbk),nb);
output_gen=zeros(NTb(tbk),ng+nrg);
for js=1:ng+nrg
    output_gen(:,js)=result.x(ins+js:Stepgd:ins+NTb(tbk)*Stepgd);   % extract the main operation data of time block tbk in the scenario node sms 
end

for jb=1:nb
 currt_dem(:,jb)=result.x(ins+in_curt+jb:Stepgd:ins+NTb(tbk)*Stepgd);
end



timeu=[staut(tbk):staut(tbk):staut(tbk)*NTb(tbk)];

figure(1)
% plot(timeu, output_gen.*100)
area(timeu,output_gen.*100)
hold on
plot(timeu,100*sum(dtn(Tsumbk+1:Tsumbk+NTb(tbk),:),2),'k','LineWidth',2)
plot(timeu,100*sum(currt_dem,2),'k--','LineWidth',2)
grid on
xlabel('time')
ylabel('MW')





