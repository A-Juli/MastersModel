%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% function to lump metabolite and compartment name
function metabolite=AddComp(met,comp)
metabolite=cell(1,36);
comp=comp(1);
for i=1:length(met)
    m=met{i};
    if ~isempty(m)
        metabolite{i}=strcat(m,'[',comp,']');
    end
end