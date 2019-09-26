%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

% fucntion to move the position of a knot in order to reduce SSR
function [step,SSRtot,SSR,to,P,co,V,C,pError]=MovePoint2(point,x,S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,max,SSR,K,order)
minim=to(point-1);
maxim=to(point+1);
step=0;
SSRtot=[x*10 SSR];
while SSRtot(end-1)>=SSR || SSR>x
    tof=to;
    tof(point)=tof(point)+(1/2)*(maxim-tof(point));
    [~,~,SSRf,~,~,~]=CalculateValues(S,R,Rr,cm,rm,t,tof,tr,stdev,stdevR,K,order);
    tor=to;
    tor(point)=tor(point)-(1/2)*(tor(point)-minim);
    [~,~,SSRr,~,~,~]=CalculateValues(S,R,Rr,cm,rm,t,tor,tr,stdev,stdevR,K,order);
    if min(SSRf,SSRr)<SSR
        if SSRf<SSR
            minim=to(point);
            to=tof;
            SSR=SSRf;
        else
            maxim=to(point);
            to=tor;
            SSR=SSRr;
        end
        SSRtot=[SSRtot SSR];     
    else
        break
    end
    if step==max
        break
    end
    step=step+1;
    fprintf('\n %1.0f\n',step);    
end
[P,co,SSR,V,C,~,pError]=CalculateValues(S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,K,order);


    