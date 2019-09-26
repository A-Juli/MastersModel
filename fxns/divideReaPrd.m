%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% to divide reaction between substrates and products
function [reactants,products]=divideReaPrd(rxn)
space = isspace(rxn);
place=findstr(1,space);
delete=place;
if ~isempty(place)
    for j=1:length(place)
        rxn(delete(j))=[];
        delete=delete-1;
    end
end
igual=findstr('=',rxn);
reactants=rxn(1:igual-1);
products=rxn(igual+1:end);