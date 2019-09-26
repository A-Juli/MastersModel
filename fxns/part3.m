%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

function Ci=part3(tj,to2,i)
if to2(i+3)-to2(i+2)~=0
    Ci=(((tj-to2(i))*((to2(i+3)-tj)^2))/((to2(i+3)-to2(i))*(to2(i+3)-to2(i+1))*(to2(i+3)-to2(i+2))))+(((to2(i+4)-tj)/(to2(i+4)-to2(i+1)))*((((tj-to2(i+1))*(to2(i+3)-tj))/((to2(i+3)-to2(i+1))*(to2(i+3)-to2(i+2))))+(((to2(i+4)-tj)*(tj-to2(i+2)))/((to2(i+4)-to2(i+2))*(to2(i+3)-to2(i+2))))));
else
    Ci=0;
end