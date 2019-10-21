%Exometabolomics Unified Procedure
%Andrei Ligema 2019 DTU
clearvars -except model_iCHOv1_S_Base
clc
close all
% Filepath input, use a spreadsheet with filepath in col 1, sheet name
% in col 2 and the experiment title in col 3
[~,paths]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\Results_Dir.xlsx');
paths=rmmissing(paths);
Num_Exps=size(paths,1);

% Module activation. 1=on 0=off
metname_formalise=0; % required for integration step, presumes all input tables share format
Outlier_Removal=1;
km_thresholding=0;
km_clustering=0;
MB_splines=1;
boundary_integration=0;

% Results Loading parameters
% start and end indices of useable data
t1=4; t2=16;
% Handle biomass as a concentration (VCD) rather than grams DW
biomass_conc=0;

% Outlier removal parameters
% Lower and upper bounds for percentile removal of points
% 0<pt_low<pt_high<100
pt_low=15; pt_high=85;

% k-means thresholding parameters

% Strictness of error convergence, 0.80<E_th<0.99
% May present as an increasing vector to evaluate multiple thresholds
E_th=0.90;
% E_th=[0.99 0.95 0.90 0.85 0.80];  
% Number of repeats for establishing 95% CI of target, more gives a better
% estimate of k for clustering
S=30;

% k-means clustering parameters
% Number of centroids, if km_thresholding is active it will override this
k_def=9;

% MB_splines parameters
HPLC=1; % is HPLC data for glucose and lactate available?
% Boundary integration parameters

% End of Parameters

% Load data
for i=1:Num_Exps
    Headers{i,1}=char(strcat(paths(i,3),"_",paths(i,2)));
%     Variables{i}='struct';
end

if metname_formalise == 1 && boundary_integration == 1
model=model_iCHOv1_S_Base;
end
    Results=table;
for I=1:Num_Exps
    clear Output Cexp Cexpl
    close all
    k=k_def;    
    Output=struct;
    spreadsheet=paths{I,1};    page=paths{I,2};
    [Cexp , Cexpl, time, Vars]=exomets_data_import (spreadsheet,page,t1,t2,biomass_conc);
    Output.CexpRaw=Cexp;    Output.Cexp=Cexp;    Output.Cexpl=Cexpl;    Output.time=time;    Output.Vars=Vars;
    
    if metname_formalise == 1 && I == 1
        % runs once per cycle to establish a table of metabolites and
        % updates the model to include any misssing metabolites
     [EMB,MetNum,model]=exomets_formal_names (model,Cexpl);
    end
    
    
    if Outlier_Removal == 1
        
    [Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pt_low,pt_high);

        Output.CexpClean=Clean_Met;
        Output.Cexp=Clean_Met;
        Vars(:,4)=Clean_Met(:,36);
        Cexp=Clean_Met;
    end
    
    if km_thresholding == 1
    [k_min,E_th]=exomets_kmeans_thresholding(E_th,S,Cexp , Cexpl, time, Vars);   
    k=k_min;
    end
    Output.K=k;
    
    if km_clustering == 1
    [Q_Pred,idx,Ret_C]=exomets_kmeans(k,Cexp , Cexpl, time, Vars, spreadsheet,page);
    Output.Q=Q_Pred;
    Output.idx=idx;
    Output.Ret_C=Ret_C;
    end
    
  
    if MB_splines == 1
     [to,P,V2,C2,SSR,MinCI2,MaxCI2,metabM,namesR,Measurements]=BDMFA_program(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2);
      Output.Measurements=Measurements; Output.to=to; Output.P=P; Output.V2=V2; Output.C2=C2; Output.SSR=SSR; Output.MinCI2=MinCI2; Output.MaxCI2=MaxCI2; Output.metabM=metabM; Output.namesR=namesR;
      time_exp=(time.*10)-(time(1)*10);
      h=figure(2);
      for J=1:length(metabM)
          subplot(6,7,J)
         plot(C2(:,J))
         hold on
         plot(time_exp(:),Measurements(:,J),'-*r')
         title(metabM(J))
      end
      figtitle=strcat(matlab.lang.makeValidName(Headers{I,1}),'_2.fig');
      savefig(h,figtitle);
%       close(h)
    end
    
    if boundary_integration == 1 && metname_formalise == 0
       
    else
        if boundary_integration == 1
        fprintf('Unable to perform boundary integration without a formal metabolite list')
        end            
    end        
  
    Results{1,I}=Output;
    Results.Properties.VariableNames{I}=matlab.lang.makeValidName(Headers{I,1});
    
    
end
    
    
