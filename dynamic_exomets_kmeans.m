clear
close all
clc

% Evaluation of steady state regions of exometabolic concentrations by
% k-mean grouping, with evaluation of fit by CHI^2
% Currently quite messy, needs cleaned up and all variables are left 
% in place for debugging
% Andrei Ligema
% 2019 DTU



[num txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_S.xlsx','Data1');
colour='b';




p=5e-4;
t1=4; t2=20;
time=num(t1:t2,1); % hours
VCD=(num(t1:t2,2));  % 10^6 cells per ml
Vol=num(t1:t2,16); % ml

xq=linspace(time(1),time(end),500);

Vars(:,1)=time;
Vars(:,2)=VCD;
Vars(:,3)=Vol;


% Importing data from Spreadsheets

%Lab daily analysis
a=[1,1;23,15];
%HPLC SUG/OA
b=[27,1;49,9];
%HPLC AA
c=[53,1;75,20];

Datag=num(a(1,1):a(2,1),6);
Addg=num(a(1,1):a(2,1),14); %after sampling
Volumeg=num(a(1,1):a(2,1),16);
GlutC(1)=Datag(1);
for i=2:length(Datag);
    GlutC(i,1)=Datag(i)-sum(Addg(1:i-1))*200/Volumeg(i);
end
GlutC=GlutC;

%HPLC
Datag=num(c(1,1):c(2,1),7);
Addg=num(a(1,1):a(2,1),14); %after sampling
Volumeg=num(a(1,1):a(2,1),16);
GlutC2(1)=Datag(1);
for i=2:length(Datag);
    GlutC2(i,1)=Datag(i)-sum(Addg(1:i-1))*200/Volumeg(i);
end
GlutC2=GlutC2;

Cexp(:,1:5)=num(t1:t2,4:8);
Cexp(:,6:14)=num(26+t1:26+t2,2:10);
Cexp(:,15:34)=num(52+t1:52+t2,2:21);

Cexpl(1,1:5)=txt(7,4+1:8+1);
Cexpl(1,6:14)=txt(33,2+1:10+1);
Cexpl(1,15:34)=txt(59,2+1:21+1);

Cexp(:,3)=GlutC(t1:t2,1)';
Cexp(:,20)=GlutC2(t1:t2,1)';


% adapted from Ivan's q calculations code, returns a list of metabolite
% concentrations (Cexp) with associated list of metabolites (Cexpl).

% Start of fitting functionality

% Create interpolated table of reactor variables





Interp_Vars(:,1)=xq;

for j=2:size(Vars,2)
    splinefun=csaps(time,Vars(:,j),p);
    for i=1:length(xq)
        Interp_Vars(i,j)=ppval(splinefun,xq(i));
    end
end
Interp_Vars(:,4)=Interp_Vars(:,2).*Interp_Vars(:,3).*0.36;
% placeholder, need to replace Col 3 here with a stepped function to
% properly reflect the sample volume 
% Gives values of time, VCD, vol and biomass for the linspace xq

for j=1:size(Interp_Vars,2)
    for i=1:length(xq)-1
        Delta_Vars(i,j)=Interp_Vars(i+1,j)-Interp_Vars(i,j);
    end
end
% produces a table of delta values for time, VCD, vol and biomass
% done establishing bioreactor variables


for j=1:size(Cexp,2)
splinefun=csaps(time,Cexp(:,j),p);
for i=1:length(xq)
    Interp_Cexp(i,j)=ppval(splinefun,xq(i));
end
end
% generates a full table of interpolated concentration values for each
% metabolite



for j=1:size(Interp_Cexp,2)
    for i=1:length(xq)-1
        DeltaC(i,j)=Interp_Cexp(i+1,j)-Interp_Cexp(i,j);
    end
end
% full table of differences for the metabolites


for j=1:size(DeltaC,2)
    for i=1:size(DeltaC,1)
        q(i,j)=DeltaC(i,j)/Interp_Vars(i,4);
    end
end
% evaluate q for the table of differences
k=1;
CHI_Flag=0;
while CHI_Flag==0

[idx,C,sumd,D]=kmeans(q,k);

for j=1:size(C,2)
    for i=1:size(idx,1)
Q_Pred(i,j)=C(idx(i),j);
    end
end
for j=1:size(q,2)
    for i=1:size(q,1)
        obs=abs(q(i,j));
        exp=abs(Q_Pred(i,j));
        Diffs(i,j)=((obs-exp)^2)/exp;
    end
end

CHI2=nansum(Diffs,'all');
dof=499-k;
IGQ=gammaincinv(CHI2,dof);



if CHI2<=0.05
    CHI_Flag=1;
else
    k=k+1;
end
end

for j=1:size(C,2)
    for i=1:size(idx,1)
Ret_DeltaC(i,j)=Q_Pred(i,j)*Interp_Vars(i,4);
    end
end

Ret_C(1,:)=Interp_Cexp(1,:);
for j=1:size(C,2)
    for i=2:size(xq,2)
        Ret_C(i,j)=Ret_C(i-1,j)+Ret_DeltaC(i-1,j);
    end
end
fprintf('Fitting complete on iteration %u \n',k);
fprintf('Quality of fit, CHI-squared p-value: %f \n', CHI2);

output_q=input('Would you like to view the fitted curves? (Y/N) \n','s');
if output_q=='Y'
    figure(1)
plot(xq(1:499),idx','o')
ylabel('Group Assignment');
xlabel('Time Interval');
ylim([0,k]);
    figure(4)
for i=1:size(Cexp,2)
    subplot(6,6,i)
    plot(xq,Interp_Cexp(:,i),'r')
    hold on
    plot(xq,Ret_C(:,i),'b')
    title(Cexpl(i))
end
else
end
clear output_q

output_q=input('Would you like to save the predicted concentration matrices? \n To save the full set of matrices for debugging / validation input F \n (Y/N/F) \n ','s');
if output_q=='Y'
    save('Predicted_Concentrations','Ret_C');
else
    if output_q=='F'
        save('Predicted_Concentractions_Debug','DeltaC','Diffs','idx','Delta_Vars','Interp_Cexp','Interp_Vars','k','Q_Pred','Ret_DeltaC','Ret_C');
    else
    end
end

output_cleanup=input('Would you like to keep the workspace after finishing? (Y/N) \n','s');
if output_cleanup~='Y'
    clear
else
end

    

% outputs:
% idx is the cluster assignment for each time span
% C is the q value for each metabolite for each of the clusters
% sumd is sum of distances to each cluster, within cluster
% D is the distances to each centroid from each point

% retreiving the metabolite concentrations from the kmeans data
% preliminary implementation, depends on all groupings being continuous

% q values are assigned to each interval


% fully retrieved set of concentrations



% validation ouputs

% figure(1)
% plot(xq(1:499),idx','o')
% ylabel('Group Assignment');
% xlabel('Time Interval');
% ylim([0,k]);

% XQ=linspace(1,500,500);
% figure(2)
% subplot(1,2,1)
% plot(xq,Interp_Cexp(:,2),'r')
% hold on
% plot(xq,Ret_C(:,2),'b')
% legend('b-spline','k-mean derived')
% xlabel('time: hrs')
% ylabel('concentration')
% title('spline vs k-mean derived concentrations')
% % ylim([0,30])
% subplot(1,2,2)
% plot(time,Cexp(:,2),'-*')
% xlabel('time: hrs')
% ylabel('concentration')
% title('raw data for concentration')
% ylim([0,30])

% figure(4)
% for i=1:size(Cexp,2)
%     subplot(6,6,i)
%     plot(xq,Interp_Cexp(:,i),'r')
%     hold on
%     plot(xq,Ret_C(:,i),'b')
%     title(Cexpl(i))
% end

% figure(3)
% CT=C';
% bar(CT)
% xlabel('Metabolite no')
% ylabel('Predicted q from k-mean')
% legend('k1','k2','k3','k4')

