%author Veronica Martinez, The University of Queensland, AIBN
%18-August-2014

% function to generate stoichiometric matrix and list of reactants
function [names,M]=CreateM(r,p,sr,sp)
names=[r(:);p(:)];
%names=DeleteEmpty(names);
names = names(~cellfun('isempty',names));
names=unique(names);
M=zeros(length(names),size(r,1));
for i=1:size(r,1);
    rsearch=r(i,:);
    for k = length(rsearch):-1:1
        if isempty(rsearch{k})
            rsearch(k)=[];
        end
    end
    for k=1:length(rsearch)
        x = strmatch(rsearch(k),names,'exact');
        if ~isempty(x)
            M(x,i)=-sr(i,k);
        end
    end
    psearch=p(i,:);
    for k = length(psearch):-1:1
        if isempty(psearch{k})
            psearch(k)=[];
        end
    end
    for k=1:length(psearch)
        x = strmatch(psearch(k),names,'exact');
        if ~isempty(x)
            M(x,i)=M(x,i)+sp(i,k);
        end
    end
end
