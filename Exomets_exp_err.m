% Experimental Error function

function [Exp_err]=Exomets_exp_err
[num , txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_S.xlsx','Data1');
Cexp1=MetsList(num,txt);
[num , txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_S.xlsx','Data2');
Cexp2=MetsList(num,txt);
[num , txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_ZeLa.xlsx','Data1');
Cexp3=MetsList(num,txt);
[num , txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_ZeLa.xlsx','Data1');
Cexp4=MetsList(num,txt);

% done compiling a list of metabolites and bioreactor variables
% Length scales
P_1=size(Cexp1,1);
P_2=P_1*2;
for i=1:size(Cexp1,2)
    clear Meas SD_Meas SD_SD SE_SD
    Meas(1:P_1,1)=Cexp1(:,i);
Meas(P_1+1:P_2,1)=Cexp3(:,i);
Meas(1:P_1,2)=Cexp2(:,i);
Meas(P_1+1:P_2,2)=Cexp4(:,i);

Meas=rmmissing(Meas);
Meas=Meas';

SD_Meas=std(Meas);
SD_SD=std(SD_Meas);
SE_SD=(SD_SD /(sqrt(length(SD_Meas))));

Exp_err(i,1)=mean(SD_Meas);
Exp_err(i,2)=SE_SD;
end



end




function Cexp=MetsList(num,txt)
t1=4; t2=20;
time=num(t1:t2,1); % hours
VCD=(num(t1:t2,2));  % 10^6 cells per ml
Vol=num(t1:t2,16); % ml
OUR=num(t1:t2,12); %mmols per L per hr

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
Cexp(:,35)=OUR;
Cexpl(1,35)={'OUR'};
Vars(:,4)=Vars(:,2).*Vars(:,3).*0.36;
Cexp(:,36)=Vars(:,4);
Cexpl(1,36)={'Biomass'};
end