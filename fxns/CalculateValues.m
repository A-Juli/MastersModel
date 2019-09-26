%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

% function to estimate mass balanced metabolic fluxes when the knots position 
% are known (to)
function [P,co,SSR,V,C,p,pError]=CalculateValues(S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,K,order)
[co,P,SSR,pError]=BSplineDMFA(K,S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,order);
V=[];
C=[];
for i=1:length(t)
    tj=t(i);
    [N, IN]=Nblendfxn(tj,to,order);
    v=calculatev(K,P,pError,N);
    c=calculateC(co,R,P,to,K,order,IN);
    V=[V; v'];
    C=[C; c'];
end
p=co;
for i=1:size(P,2)
    p=[p;P(:,i)];
end

