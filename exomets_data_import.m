function [Cexp , Cexpl, time, Vars]=exomets_data_import(spreadsheet,page,t1,t2)

% [num , txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_S.xlsx','Data1');
[num , txt]=xlsread(spreadsheet,page);



% t1=4; t2=20;
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
% adapted from Ivan's q calculations code, returns a list of metabolite
% concentrations (Cexp) with associated list of metabolites (Cexpl).
