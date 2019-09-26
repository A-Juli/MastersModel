%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

function Ai=part1(tj,to,i)
Ai=((tj-to(i))^3)/((to(i+3)-to(i))*(to(i+2)-to(i))*(to(i+1)-to(i)));
