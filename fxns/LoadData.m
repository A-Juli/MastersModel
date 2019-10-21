%author Veronica Martinez, The University of Queensland, AIBN
%13-March-2014
%modified version by Andrei Ligema September 2019 DTU Biosustain

% function to load experimental data and a metabolic model

function [t,tr,cm,rm,stdev,stdevR,R,Rr,S,metabM,namesR,Measurements]=LoadData(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2,biomass)

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
Measurements(:,2)=Vars(:,2); % VCD
Measurements(:,3)=Cexp(:,5); % ammonium
Measurements(:,4)=Cexp(:,36); % Biomass
Measurements(:,5:31)=Cexp(:,8:34); % Pyruvate to Proline
if HPLC == 0 % Glucose and Lactate
Measurements(:,32)=Cexp(:,1); 
Measurements(:,33)=Cexp(:,2); 
else
Measurements(:,32)=Cexp(:,6); 
Measurements(:,33)=Cexp(:,7); 
end
Measurements(:,34:35)=Cexp(:,37:38); % Na and K
[~ , txt]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\IMM_q_calculation_yields_DMFA_S.xlsx','Met_List');
metabM=txt'; % time



% Interpolating through missing data points
for j=1:size(metabM,2)
   
    interp_met=Measurements(:,j); % import raw/filtered data
    missing=find(isnan(interp_met)); % check for absent datapoints
    if isempty(missing)==0         
    t_missing=time(missing);      % index missing datapoints by time
    interp_met(:,2)=time;         % add time index to points  
    interp_met(:,3)=[1:1:size(interp_met,1)]; % add global index to points
     % add space for flags
    interp_met=rmmissing(interp_met); % cull NaN points (listed in t_missing by time)
    Exit_flag=0;                  % exit state for the while loop
    k=1; % position in metabolite 
    x=3; % interpolation segment
while Exit_flag==0
    if interp_met(k+1,3)~=(interp_met(k,3)+1) % check for adjacency of points
             x=x+1;                           % create a new interpolation segment
             interp_met(1:size(interp_met,1),x)=zeros;
             interp_met(k,x)=1; interp_met(k+1,x)=1; % flag the start and end of the segment
    end
    k=k+1;   % advance to the next point
    if k==size(interp_met,1) % reached the end of the array
        Exit_flag=1;
    end
end

if x~=3 % if at least one gap was detected
for y=4:x 
    clear interp_stage interp_local_time interp_local_vals
    ITP=find(interp_met(:,y)); % The points associated with this interpolation segment by local indexing
    interp_local_time=interp_met(ITP(1):ITP(2),2);
    interp_local_vals=interp_met(ITP(1):ITP(2),1);
    interp_fn=griddedInterpolant(interp_local_time,interp_local_vals); % create interp fn for this segment
    T1=interp_met(ITP(1),3); T2=interp_met(ITP(2),3); % identify the start and end points by global indexing
    T_span=(T2-T1)-1; 
    T_interp=t_missing(1:T_span);
    t_missing(1:T_span)=[];
    
    interp_stage=interp_fn(T_interp); % recover interpolated values from the function
    if y==4
        interp_vals=interp_stage;
    else
        interp_vals=cat(1,interp_vals,interp_stage);
    end
    
end
Measurements(missing,j)=interp_vals;

end    
    
    
    end
     clear interp_met missing t_missing interp_fn interp_vals
end     
 
% removing unused metabolites, temporary
Rem=[5 6 7 8 9 10 11 14 34 35];
metabM(Rem)=[];
Measurements(:,Rem)=[];

t=Measurements(:,1);
Measurements(:,1:2)=[];
cm=Measurements;
metabM(1:2)=[];

stdev=cm*0.05;% 5% error of measurements
stdev(find(stdev<0.01))=0.01;
% rates data
rm=(Cexp(:,35).*(Vars(:,3))./1000);

for j=1:size(rm,2)
    interp_met=rm(:,j);
    missing=find(isnan(interp_met));
    if isempty(missing)==0 
    t_missing=time(missing);
    interp_met(:,2)=time;
    interp_met=rmmissing(interp_met);
    interp_fn=griddedInterpolant(interp_met(:,2),interp_met(:,1));
    interp_vals=interp_fn(t_missing);
    rm(missing,j)=interp_vals;
    end
     clear interp_met missing t_missing interp_fn interp_vals
end     
 

stdevR=rm*0.05;;
stdevR(find(stdevR<0.01))=0.01;

namesR={'Oxygen[e]'};
tr=t;

%Load metabolic model
[~, Reactions]=xlsread('C:\Users\anjuli\Documents\B-DMFA - generic order\TemperatureShift.xlsx','Model');
[~,Mc,~,~,~,namesc,transport]=reactions(Reactions);
R=[];
for i=1:length(metabM)
    prePosition=startsWith(namesc,metabM{i});
    prePosition=find(prePosition);
    R(i,:)=Mc(prePosition,:);
end
for i=1:length(namesR)
    prePosition=startsWith(namesc,namesR{i});
    prePosition=find(prePosition);
    Rr(i,:)=Mc(prePosition,:);
end


%S stochiomatric matrix of balanced metabolites
[S,Names]=DeleteRows(namesc,Mc,'e');
[S,Names]=DeleteRows(Names,S,biomass);

%check data in figure 1
figure(1)
for i=1:size(cm,2);
    subplot(6,6,i),plot(t,cm(:,i));title(metabM{i});
end

