%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

function [co,P,SSR,pError]=BSplineDMFA(K,S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,order)
% Dynamic metabolic flux analysis (DMFA) using a B-spline fitting of any
% order
Me=IntN(to,order);
H=zeros(size(S(1,:)*K,2)*(length(to)-order)+size(cm,2)); 
J=zeros(size(S(1,:)*K,2)*(length(to)-order)+size(cm,2),1);
for j=1:length(t)
    tj=t(j);
    [N, IN]=Nblendfxn(tj,to,order);
    IN=Me*IN;
    for i=1:size(R,1)
        Si=R(i,:);
        diagonal=Si*K;
        Di=zeros(length(to)-order,(length(to)-order)*length(diagonal));
        for l=1:size(Di,1)
            Di(l,(l-1)*length(diagonal)+1:l*length(diagonal))=diagonal;
        end
        Ei=zeros(1,size(cm,2));
        Ei(i)=1;
        dcidp=[Ei (IN'*Di)];
        Wi=1/(stdev(j,i)^2);
        H=H+dcidp'*Wi*dcidp;
        J=J+dcidp'*Wi*cm(j,i);
    end
    if ~isempty(find(tr==tj, 1))
        for i=1:size(Rr,1)
            Si=Rr(i,:);
            diagonal=Si*K;
            Di=zeros(length(to)-order,(length(to)-order)*length(diagonal));
            for l=1:size(Di,1)
                Di(l,(l-1)*length(diagonal)+1:l*length(diagonal))=diagonal;
            end
            dridp=[zeros(1,size(R,1)) (N'*Di)];
            Wir=1/(stdevR(find(tr==tj),i)^2);
            H=H+dridp'*Wir*dridp;
            J=J+dridp'*Wir*rm(find(tr==tj),i);
        end
    end
end
H1=CalculationInverse(H);        
x=H1*J;
co=x(1:size(cm,2));
p=x(size(cm,2)+1:end);
P=[];
dim=size(K,2);
for i=1:length(to)-order
    P=[P p((i-1)*dim+1:i*dim)];
end
SSR=0;
CovCoP=H1;
CovP=CovCoP(size(co,1)+1:end,size(co,1)+1:end);
dim=size(P,2);
dim2=size(P,1);
pError=[];
for i=1:dim
    j=i*dim2-(dim2-1);
    Covp=CovP(j:j+dim2-1,j:j+dim2-1);
    pError(:,i)=sqrt(diag(K*Covp*K'));
end
for j=1:length(t)
    tj=t(j);
    [N, IN]=Nblendfxn(tj,to,order);
    for i=1:size(cm,2)
        Si=R(i,:);
        ci=calculateC(co(i),Si,P,to,K,order,IN);
        SSR=SSR+((ci-cm(j,i))^2)/(stdev(j,i)^2);
    end
    for i=1:size(rm,2)
        ri=Rr(i,:)*calculatev(K,P,pError,N);
        SSR=SSR+((ri-rm(find(tr==tj),i))^2)/(stdevR(find(tr==tj),i)^2);
    end
end










