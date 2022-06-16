


% computing the the sets of all parents nodes of m  until the stage
    
function [p0m,pim]=phim(NM,em,fathernode,Nodeprob)

% p0m{m}       -    It contains the set of all parents node of m including m iteself
% pim(m)       -    Probability that tree node m occurs



for m=1:NM
    
Epm=em(m);
pm=[];
pm=m*ones(1,Epm);
prm=zeros(1,Epm-1);
for i=1:Epm-1
  pm(1,Epm-i)=fathernode(1,pm(1,Epm-i+1));
  prm(1,Epm-i)=Nodeprob(pm(1,Epm-i),pm(1,Epm-i+1));
end

p0m{m}=pm;

%Probability that tree node m occurs

pim(m)=prod(prm);

end