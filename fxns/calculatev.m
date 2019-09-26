%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

%function to estimate fluxes knowing the control points (P) and knots (to)
function [v,Stdv]=calculatev(K,P,preError,N)
v=K*P*N;
Stdv=preError*(N.^2);
