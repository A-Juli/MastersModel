%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

function Di=part4(tj,to2,i)
if to2(i+4)-to2(i+1)~=0 && to2(i+4)-to2(i+3)~=0 && to2(i+4)-to2(i+3)~=0
    Di=((to2(i+4)-tj)^3)/((to2(i+4)-to2(i+1))*(to2(i+4)-to2(i+2))*(to2(i+4)-to2(i+3)));
else
    Di=0;
end