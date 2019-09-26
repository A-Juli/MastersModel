%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% to divide substrates or products part of the reaction between 
% metabolites and stoichiometric coeficients
function [met,smet]=divideMetab(part)
plus=findstr('+',part);
l=length(plus)+1;
smet=zeros(1,36);
met=cell(1,36);
 for j=1:l
      if j==1
         if isempty(plus)
           m=part(1:end);
         else
             m=part(1:plus(j)-1);
         end
       else 
         if j<l
            m=part(plus(j-1)+1:plus(j)-1);
         else
             m=part(plus(j-1)+1:end);
         end
      end
      smet(1,j)=1;
      pa1=findstr('(',m);
      pa2=findstr(')',m);
       if ~isempty(pa1) && pa1(1)==1
           smet(1,j)=str2double(m(2:pa2(1)-1));
           m=m(pa2(1)+1:end);
       end 
       met{1,j}= m;
 end