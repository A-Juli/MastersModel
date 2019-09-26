%Exometabolomics Unified Procedure
%Andrei Ligema 2019 DTU
clear 
clc
close all
% Filepath input, use a spreadsheer with filepath in col 1 and sheet name
% in col 2
[~,paths]=xlsread('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\Results_Dir.xlsx');
Num_Exps=size(paths,1);

% Module activation. 1=on 0=off
Outlier_Removal=1;
km_thresholding=0;
km_clustering=0;
MB_splines=1;
boundary_integration=0;

% Results Loading parameters
% start and end indices of useable data
t1=4; t2=20;

% Outlier removal parameters
% Lower and upper bounds for percentile removal of points
% 0<pt_low<pt_high<100
pt_low=25; pt_high=75;

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

% Boundary integration parameters

% End of Parameters

% Load data
for i=1:Num_Exps
    Headers{i,1}=char(strcat(paths(i,3),"_",paths(i,2)));
%     Variables{i}='struct';
end
    Results=table;
for I=1:Num_Exps
    clear Output
    k=k_def;    
    Output=struct;
    spreadsheet=paths{I,1};    page=paths{I,2};
    [Cexp , Cexpl, time, Vars]=exomets_data_import(spreadsheet,page,t1,t2);
    Output.CexpRaw=Cexp;    Output.Cexp=Cexp;    Output.Cexpl=Cexpl;    Output.time=time;    Output.Vars=Vars;
    if Outlier_Removal == 1
        
    [Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pt_low,pt_high);
        % add linear interpolation to reconstruct a Cexp matrix here? 
        % only needed for input to veronica's code and common to basic Cexp
        % so just run at BDMFA step
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
     [to,P,V2,C2,SSR,MinCI2,MaxCI2,metabM,namesR]=BDMFA_program(Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2);
      Output.to=to; Output.P=P; Output.V2=V2; Output.C2=C2; Output.SSR=SSR; Output.MinCI2=MinCI2; Output.MaxCI2=MaxCI2; Output.metabM=metabM; Output.namesR=namesR;
      
    end
    if boundary_integration == 1
        % rework Exomets Integration to unify with Cexp format
    end        
  
    Results{1,I}=Output;
    Results.Properties.VariableNames{I}=Headers{I,1};
end
    
    
