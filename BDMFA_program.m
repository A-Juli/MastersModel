function [to,P,V2,C2,SSR,MinCI2,MaxCI2,metabM,namesR,Measurements]=BDMFA_program(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2)

%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

%B-DMFA program to estimate mass balance fluxes, the main outputs are: 
% to, the knots positions
% P, control points
% V2, estimated fluxes over time
% C2, fitted concentrations over time
% SSR, the variance-weighted sum of squared residuals 
% MinCI2, min rates CI
% MaxCI2, max rates CI 

% clear ;clc;
tic
% load data cm ( measured concentrations), rm (measured rates),
% S (matrix with balanced metabolites), R (concentrations measurement matrix), 
% t (time points of measured data), stdev (standar desviation of data, concentrations)
% Rr (rates measurement matrix), stdevR (standar desviation of data, rates)
addpath ./fxns 

[t,tr,cm,rm,stdev,stdevR,R,Rr,S,metabM,namesR,Measurements]=LoadData(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2,'b');
% run YeastExample 

order=3;
to=[t(1) t(end)]';
to= [ones(1,order-1)*min(to) to' ones(1,order-1)*max(to)]; % with (order-1) extra knots at the extreme time points
K=null(S);
[P,co,SSR,~,~,p,pError]=CalculateValues(S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,K,order);
SSRi=SSR;

% Start with the algoritm to locate the internal knots
% Estimate goodness-of-fit, with chi-squared considering 95% CI
n= length(t)*(size(cm,2))+length(tr)*(size(rm,2)); % fitted measurements
pa=length(p); % estimated parameters
np=n-pa; % redundant data
x=chi2inv(1-0.025,np); % 95% CI

Step=0;
divisions=1;
while SSR>x && np > size(P,1)+1
    toint=to(order:end-order+1);
    SSRcheck=zeros(divisions,1);
    for i=1:divisions
        if i==1
            mini=1;
            maxi=max(find(t<=toint(2)));
        elseif i==length(divisions)
            mini=min(find(t>=toint(end-1)));
            maxi=legth(t);
        else
            mini=min(find(t>=toint(i)));
            maxi=max(find(t<=toint(i+1)));
        end
        SSRcheck(i)=EstimateSSR(to,t,tr,cm,rm,R,Rr,P,K,co,stdev,stdevR,mini,maxi,pError,order);
    end
    workplace=find(SSRcheck==max(SSRcheck));
    extraP=toint(workplace)+(1/2)*(toint(workplace+1)-toint(workplace));
    to=[to(1:workplace+order-1) extraP to(workplace+order:end)];
    [P,~,SSR,~,~,p,pError]=CalculateValues(S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,K,order);
    pa=length(p)+length(to)-6; % estimated parameters
    np=n-pa; % redundant data
    x=chi2inv(1-0.025,np); %95% CI
    % move the position of extra knot
    [step,SSRtot,SSR,to,P,co,~,~,pError]=MovePoint2(find(to==extraP),x,S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,100,SSR,K,order);
    %move all internal knots 
    if Step>0
        SSRpre=SSR+1;
        while SSR<SSRpre
            SSRpre=SSR;
            for i=1:length(to)-(order)*2
                [~,~,SSR,to,P,co,~,~,pError]=MovePoint2(i+order,x,S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,100,SSR,K,order);
            end
        end
    end
    Step=Step+1;
    divisions=divisions+1;
    fprintf('\n %1.0f\n',Step);
end

%The final result (fluxes)
[P,co,SSR,V,C,p,pError]=CalculateValues(S,R,Rr,cm,rm,t,to,tr,stdev,stdevR,K,order);
Error=[];
for i=1:length(t)
    ti=t(i);
    [N, IN]=Nblendfxn(ti,to,order);
    [v,Stv]=calculatev(K,P,pError,N);
    Error=[Error; Stv'];
end
% determine 95% confidence interval
MinCI=V-(2*Error);
MaxCI=V+(2*Error);

t2=t(1):0.1:t(end);
V2=[];
C2=[];
Error2=[];
Me=IntN(to,order);
for j=1:length(t2)
    tj=t2(j);
    [N, IN]=Nblendfxn(tj,to,order);
    v=K*P*N;
    c=co+(R*K*P*Me*IN);
    N=N';
    Stdv2=pError*(N.^2)';
    V2=[V2; v'];
    C2=[C2; c'];
    Error2=[Error2; Stdv2'];
end
% determine 95% CI
MinCI2=V2-(2*Error2);
MaxCI2=V2+(2*Error2);

toc
% run PlotYeast

