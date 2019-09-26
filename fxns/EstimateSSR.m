%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

%function to estimate the variance-weighted sum of squared residuals (SSR)
%knowing the fitting curve
function SSR=EstimateSSR(to,t,tr,cm,rm,R,Rr,P,K,co,stdev,stdevR,mini,maxi,preError,order)
SSR=0;
SSR2=[];
for j=mini:maxi
    tj=t(j);
    [N, IN]=Nblendfxn(tj,to,order);
    for i=1:size(cm,2)
        Si=R(i,:);
        ci=calculateC(co(i),Si,P,to,K,order,IN);
        SSR=SSR+((ci-cm(j,i))^2)/(stdev(j,i)^2);
    end
    if ~isempty(find(tr==tj, 1))
        for i=1:size(rm,2)
            ri=Rr(i,:)*calculatev(K,P,preError,N);
            SSR=SSR+((ri-rm(find(tr==tj),i))^2)/(stdevR(find(tr==tj),i)^2);
        end
    end
    SSR2(j)=SSR;
end
SSR=sum(SSR2);