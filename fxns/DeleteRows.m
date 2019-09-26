%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% to remove boundary metabolites from stoichiometric matrix
function [Mc_f,namesc_f]=DeleteRows(namesc_b,Mc_b,b)
delete=[];
for  i = 1:length(namesc_b)
    met=namesc_b{i};
    starting=findstr('[',met);
    comp=met(starting+1:end-1);
    flag=strmatch(b,comp,'exact');
    if ~isempty(flag)
        delete=[delete;i];
    end
end
namesc_f=namesc_b;
namesc_f(delete)=[];
Mc_f=Mc_b;
Mc_f(delete,:)=[];
