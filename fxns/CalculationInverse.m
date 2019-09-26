%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

function H1=CalculationInverse(H)
if cond(H)>1e10
    [U,S2,V]=svd(H);
    S2(find(S2<1e-13))=0;
    check=diag(S2);
    check(find(check==0))=[];
    check=1./check ;
    for i=1:length(check)
        S2(i,i)=check(i);
    end
    S2=transpose(S2);
    H1=V*S2*(U^-1);
else
    H1=H^-1;
end