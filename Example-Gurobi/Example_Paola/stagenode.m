function  em=stagenode(NE,Nem)


% em - the entry i contains the epoch to which the node i belongs to

em=[];
for i=1:NE

    em=[em,i*ones(1,Nem(1,i))];

end