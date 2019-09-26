%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% to remove compartment name of metabolite name
function [metabolite,comp]=RemoveComp(met)
metabolite=cell(1,36);
comp=cell(1,36);
for i=1:length(met)
    m=met{i};
    if ~isempty(m)
        comp{1,i}=m(end-1:end-1);
        metabolite{i}=m(1:end-3);
    end
end