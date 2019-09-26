%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2014

%function to estimate the fitting concentrations, knowing control points(P) 
%and knots (to)
function c=calculateC(co,R,P,to,K,order,IN)
Me=IntN(to,order);
c=co+(R*K*P*Me*IN);
    

    