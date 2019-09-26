%author Veronica Martinez, The University of Queensland, AIBN
%13-March-2014
%modified version by Andrei Ligema September 2019 DTU Biosustain

% function to load experimental data and a metabolic model

function [t,tr,cm,rm,stdev,stdevR,R,Rr,S,metabM,namesR]=LoadData(Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2,biomass)

% sheet='RawData'; %strcat('RawData',num2str(x));
% [num , txt]= xlsread(excelfile,sheet);
% % [Measurements metabM]= xlsread(excelfile,sheet);
% t1=4; t2=16;
% % Temporary implementation, this may be replaced with a trimmed data set
% % based on outlier removal 
% metabM(1)=txt(2,16);
% metabM(2)=txt(2,17);
% metabM(3)=txt(2,41);
% metabM(4)={'Biomass[b]'};
% for i=1:22
%     metabM(i+4)=txt(2,(61+(4*(i-1))));
% end
% 
% 
% Measurements(:,1)=num(t1:t2,16);
% Measurements(:,2)=num(t1:t2,19);
% Measurements(:,3)=num(t1:t2,43);
% Measurements(:,4)=Measurements(:,2).*num(t1:t2,4).*0.36;
% for i=1:22
%     Measurements(:,i+4)=num(t1:t2,(63+(4*(i-1))));
% end

% finished importing data from Ivan's sheets

% New implementation, import data from previous modules
Measurements(:,1)=time;
Measurements(:,2)=Vars(:,2);
Measurements(:,3)=Cexp(:,5);
Measurements(:,4)=Cexp(:,36);
Measurements(:,5:24)=Cexp(:,15:34);
Measurements(:,25)=Cexp(:,6);
Measurements(:,26)=Cexp(:,7);

[~ , txt]=xlsread(spreadsheet,'RawData');
metabM(1)=txt(2,16);
metabM(2)=txt(2,17);
metabM(3)=txt(2,41);
metabM(4)={'Biomass[b]'};
for i=1:22
    metabM(i+4)=txt(2,(61+(4*(i-1))));
end




% Improved solution would involve the lookup method, need to standardise
% how metabolites are listed on spreadsheets to make this more effective

% Interpolating through missing data points
for j=1:size(metabM,2)
    interp_met=Measurements(:,j);
    missing=find(isnan(interp_met));
    if isempty(missing)==0
    t_missing=time(missing);
    interp_met(:,2)=time;
    interp_met=rmmissing(interp_met);
    interp_fn=griddedInterpolant(interp_met(:,2),interp_met(:,1));
    interp_vals=interp_fn(t_missing);
    Measurements(missing,j)=interp_vals;
    end
     clear interp_met missing t_missing interp_fn interp_vals
end     
            

t=Measurements(:,1);
Measurements(:,1:2)=[];
cm=Measurements;
metabM(1:2)=[];

stdev=cm*0.05;% 5% error of measurements
stdev(find(stdev<0.01))=0.01;
%no rates data
rm=[];stdevR=[];Rr=[];namesR=[];tr=[];

%Load metabolic model % introduce step here to amend the model for missing
%metabolites
[~, Reactions]=xlsread('C:\Users\anjuli\Documents\B-DMFA - generic order\TemperatureShift.xlsx','Model');
[~,Mc,~,~,~,namesc,transport]=reactions(Reactions);
R=[];
for i=1:length(metabM)
    prePosition=startsWith(namesc,metabM{i});
    prePosition=find(prePosition);
    R(i,:)=Mc(prePosition,:);
end

%S stochiomatric matrix of balanced metabolites
[S,Names]=DeleteRows(namesc,Mc,'e');
[S,Names]=DeleteRows(Names,S,biomass);

%check data in figure 1
figure(1)
for i=1:size(cm,2);
    subplot(6,6,i),plot(t,cm(:,i));title(metabM{i});
end

