%Exometabolomics Unified Procedure
%Andrei Ligema 2019 Centre for Biosustainability
clearvars -except model_iCHOv1_S_Base model_BDMFA_Base compSymbolList compNameList EMB MetNum model Cexp_list
clc
close all
% Filepath input, use a spreadsheet with filepath in col 1, sheet name
% in col 2, the experiment title in col 3 and the experiment index number
% in col 4
[groups,paths]=xlsread('C:\Users\Andre\Documents\MATLAB\cobratoolbox\Results_Dir.xlsx');
paths=rmmissing(paths);
Num_Exps=size(paths,1);

% setting up grouped experiments for averaging of results
Exp_unique=unique(groups);
for g=1:length(Exp_unique)
   
    g_assignment(:,g)=find(groups==Exp_unique(g));
    
end

% Module activation. 1=on 0=off
% Universal
cobra_init=0; % Initialises COBRA, required for most stages, may be run once at start of work
model_import=1; % Imports the GSM used, requires SBML in .xml format
metname_formalise=0; % required for integration & GMEA presumes all input tables share format
                     % constructs a list of metabolites fully consistent
                     % with the used GSM. May be run once as long as the
                     % output variables are not cleared between executions
% Per-repeat
Outlier_Removal=1;  % Attempts to remove outliers from isolated datasets by gradient prediction
                    % and frequency domain analysis. Currently integral to
                    % the procedure. If not desired then a function should
                    % be inserted where outlier_removal is called,
                    % replacing its outputs
% Per-experiment
km_thresholding=0;  % Aggregates biological replicates to determine the minimum number of clusters
                    % required to model the experimental data to a given
                    % degree of accuracy
km_clustering=0;    % segregates an experimental run into k clusters and generates a table of predicted 
                    % flux values for each measured metabolite within each
                    % time region
GMEA_process=0;     % Gross Measurement Error Analysis: adapted from Quek et al. 2010
% Per-repeat
MB_splines=0;       % BDMFA procedure, adapted from Veronica
boundary_integration=0; % Integrates flux predictions from km_clustering or GMEA to produce a 
                        % set of updated GSMs in SBML format with new
                        % objective and boundary conditions



% Results Loading parameters
% start and end indices of useable data
t1=4; t2=22;

% Locations in the spreadsheet for Metabolite timecourse, see exomets_data_import for guidance on this 
%Lab daily analysis
a=[1,1;23,15];
%HPLC SUG/OA
b=[27,1;49,9];
%HPLC AA
c=[53,1;(53+t2-1),20];
% Biomass index within measurements list
BM_idx=36; % not currently utilised


% Handle biomass as a concentration (VCD) rather than grams DW
biomass_conc=0;

% Name Formalisation Parameters
% External metabolite tag
ExTag='ext_b';
% If not all metabolites measured are available in the model then specify
% the used ones here
Util_mets=[5 6 7 12 15:34 36];
if isempty(Util_mets)==1
    warning('No involved metabolites specified, default is to assume ALL are present in the used model. If this is not the case please update the variable "Util_mets" with the indices of the utilised metabolites ');
end
    

% Outlier removal parameters

pre_gradient=0;
pre_gradient_range=5;
exp_weight=1;
Exp_Err=0.1;

% k-means thresholding parameters
% Strictness of error convergence, 0.80<E_th<0.99

E_th=0.95;
 
% Number of repeats for establishing 95% CI of target, more gives a better
% estimate of k for clustering
S=30;

% k-means clustering parameters
% Number of centroids, if km_thresholding is active it will override this
k_def=7;

% MB_splines parameters
HPLC=1; % is HPLC data for glucose and lactate available?
% Boundary integration parameters

% End of Parameters

% Load data
for i=1:Num_Exps
    Headers{i,1}=char(strcat(paths(i,3),"_",paths(i,2)));
%     Variables{i}='struct';
end

    if cobra_init == 1 
     initCobraToolbox   
    end
    
    if model_import == 1 
        modelDirectory='C:\Users\Andre\Documents\MATLAB\cobratoolbox\ExoMets';
        modelFileName_base='ivan.v1.xml';
        modelFileName_base = [modelDirectory filesep modelFileName_base];
        model = readCbModel(modelFileName_base);
        compSymbolList=model.comps;
        compNameList=model.compNames;
    end
    


%     Results=table; % depricated since introduction of parallel
%     replicates, I should reinstate this at some point
%     Intention was to output results of each replicate/state as a
%     structure containing all the relevant progress and final variables in
%     order to allow for use/debugging without the need to pause
%     in-progress.
for I=1:Exp_unique
    
    % replicates marked with the same experiment number are run in parallel
    % and share error estimation data
    
    
    exp_index=g_assignment(:,I); 
    
    
    clear Output Cexp Cexpl
    close all
    k=k_def; % relevant only if k_means thresholding is turned off, not relevant unless user wants a specific number of clusters
    
           
    
    % Output=struct; % depricated since introduction of parallel
    % replicates, I should reinstate this at some point
    
    
    % load experiment metadata: 
     spreadsheet=paths{exp_index(1),1};    page=paths{exp_index(1),2};
     [Cexp , Cexpl, time, Vars]=exomets_data_import (spreadsheet,page,t1,t2,biomass_conc,a,b,c);
    %     Output.CexpRaw=Cexp;    Output.Cexp=Cexp;    Output.Cexpl=Cexpl;    Output.time=time;    Output.Vars=Vars;
    

    
    
    % runs once per experiment to formalise naming conventions
    % a LOT of variable nightmares start here due to the
            % mandatory-but-not-really nature of this function
           
    if metname_formalise == 1 
    [EMB,MetNum,model,Cexp_list]=experiment_initial(metname_formalise,paths,exp_index,t1,t2,biomass_conc,model,ExTag,Util_mets,a,b,c);
    else
        % set these index variables to default state utilising all
        % metabolites
         EMB=Cexpl;
        Cexp_list=[1:length(Cexpl)];
         MetNum=[]; % this is fine being empty as long as the user doesn't attempt to interface with a model at any point
    
    end
    
    % Outlier removal
    % runs replicates individually to clear outliers   
        
        Cexp_clean=table;
        Vars_clean=table;
        for U=1:length(exp_index)  
     [Clean_Met,Vars]=outlier_cleaning(exp_index,paths,pre_gradient,pre_gradient_range,exp_weight,Exp_Err,t1,t2,biomass_conc,U,a,b,c,Outlier_Removal);
        Cexp_clean{1,U}={Clean_Met}; 
        Vars_clean{1,U}={Vars};
        
        end

   
   % clusters replicates to evaluate 
    
    
    if km_thresholding == 1 || km_clustering == 1
 
[Q_tables,Q_rates_table,LB_table,UB_table,Q_rates_GMEA,revisedQ]=experiment_clustered(km_thresholding,km_clustering,k,Cexp_clean,Vars_clean,E_th,Cexpl,exp_index,Cexp_list,S,time,model,EMB,MetNum,GMEA_process);

    end
    
 
    
    %%%%%%%%%%
  if MB_splines == 1 || boundary_integration == 1 
     
   metname_override=1;   % override switch for cases where the outputs of metname_formalise are added manually by the user or they are carried from a previous run.
      
      if metname_formalise==1 || metname_override==1
      for U=1:length(exp_index)
      Cexp=Cexp_clean{1,U}{1,1};
      Vars=Vars_clean{1,U}{1,1};
     spreadsheet=paths{exp_index(U),1};    page=paths{exp_index(U),2};
     LB=LB_table{1,U}{1,1}; UB=UB_table{1,U}{1,1};
   
      experiment_outcomes(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2,model,EMB,MetNum,LB,UB,compSymbolList,compNameList,MB_splines,boundary_integration,Headers);
     
      end
       else
          error('Cannot interact with models properly without a formal list of metabolite names and index locations, please eneable metname_formalise or its associated override if a manually inserted list is used');
      end
  end
    
    
    

  
%     Results{1,I}=Output;
%     Results.Properties.VariableNames{I}=matlab.lang.makeValidName(Headers{I,1});
%     

    
end
    
function [EMB,MetNum,model,Cexp_list]=experiment_initial(metname_formalise,paths,exp_index,t1,t2,biomass_conc,model,ExTag,Util_mets,a,b,c)
    spreadsheet=paths{exp_index(1),1};    page=paths{exp_index(1),2};
    [Cexp , Cexpl, time, Vars]=exomets_data_import (spreadsheet,page,t1,t2,biomass_conc,a,b,c);
    
        % runs once per cycle to establish a table of metabolites and
        % updates the model to include any misssing metabolites
        if isempty(Util_mets)==1
            Util_mets=[1:length(Cexpl)];
        end
     [EMB,MetNum,model,Cexp_list]=exomets_formal_names (model,Cexpl,ExTag,Util_mets);
     model=model;
    
    
end

function [Clean_Met,Vars]=outlier_cleaning(exp_index,paths,pre_gradient,pre_gradient_range,exp_weight,Exp_Err,t1,t2,biomass_conc,U,a,b,c,Outlier_Removal)
    spreadsheet=paths{exp_index(U),1};    page=paths{exp_index(U),2};
    [Cexp , Cexpl, time, Vars]=exomets_data_import (spreadsheet,page,t1,t2,biomass_conc,a,b,c);  

if Outlier_Removal==1
[Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pre_gradient,pre_gradient_range,exp_weight,Exp_Err);

        Output.CexpClean=Clean_Met;
        Output.Cexp=Clean_Met; % this replicates the cleaned metabolite data for a purpose I can't quite place
        Vars(:,4)=Clean_Met(:,36); % need to replace this indexing to use an input variable instead of relying on a specific reference to Ivan's experiment layout
      %  Clean_Met; % Not sure what this is doing or why
else
    Output.CexpClean=Cexp;
    Output.Cexp=Cexp;
end
end

function [Q_tables,Q_rates_table,LB_table,UB_table,Q_rates_GMEA,revisedQ]=experiment_clustered(km_thresholding,km_clustering,k,Cexp_clean,Vars_clean,E_th,Cexpl,exp_index,Cexp_list,S,time,model,EMB,MetNum,GMEA_process)

if km_thresholding==1
    [k_min,E_th]=exomets_kmeans_thresholding(E_th,S,Cexp_clean , Cexpl, time, Vars_clean);   
    k=k_min;  
end

if km_clustering==1
Q_tables=table;
Q_rates_table=table;
LB_table=table;
UB_table=table;
Int_Vars=table;

for B=1:length(exp_index)

    [Q_Pred,idx,Ret_C,IDX_lin_it,xq]=exomets_kmeans(k,Cexp_clean , Cexpl, time, Vars_clean,B);
IDX_lin(:,B)=IDX_lin_it;



% these tables replace the functionality of previous variables
% they are arranged in order of experimental repeats by column, simply
% extract the n'th index to retrieve the data for other modules where
% required
% Cexp_list(end+1)=36; % adding biomass to the listing


    Q_tables{1,B}={Q_Pred};
    Q_temp=Q_Pred;
    for r=1:k
       Q_mean=Q_temp(IDX_lin_it==r,:);
       Q_mean=(Q_mean(1,:)).*-1;
       Q_rates_hold(r,:)=Q_mean;       
       LB_table_hold(r,:)=Q_mean-abs(Q_mean.*0.1);       
       UB_table_hold(r,:)=Q_mean+abs(Q_mean.*0.1);
       Biom_Q(r,B)=Q_mean(1,36);
    end
    Q_rates_table{1,B}={Q_rates_hold};
    Q_rates_GMEA{1,B}={Q_rates_hold(:,Cexp_list)};
    LB_table{1,B}={LB_table_hold(:,Cexp_list)};
    UB_table{1,B}={UB_table_hold(:,Cexp_list)};
    
    
end

end % end of kmeans clustering



% generate standard error tables for experiments
for p=1:size(Cexpl,2)
    clear Meas_err Meas
    for s=1:size(Cexp_clean,2)
    Meas(:,s)=Cexp_clean{1,s}{1,1}(:,p);
    Mu(:,s)=Vars_clean{1,s}{1,1}(:,4);
    Meas(:,s)=Meas(:,s).*Mu(:,s);
    end
    Meas=rmmissing(Meas);
   for r=1:size(Meas,1)
       Meas_err(r)=(std(Meas(r,:)));
       %Meas_err(r)=abs((((Meas(r)*1.1)-Meas(r))^2)/Meas(r));      
   end
   Meas_err=rmmissing(Meas_err);
   if isempty(Meas_err)==1
      Meas_err=0; 
   end
   se(p,1)=mean(Meas_err)/sqrt(length(Meas_err));   
end

se=se(Cexp_list,1);
revisedQ=table;

if GMEA_process ==1
for ver=1:length(exp_index)

    Q_rates=Q_rates_GMEA{1,ver}{1,1};
    Biom=Biom_Q(:,ver);
for U=1:k  
    
[dataOUT,rates_hat,F_hat]=GMEA(model,Q_rates,se,MetNum,EMB,U,Biom);
revisedQ{ver,U}=dataOUT;

Q_bal(U,:)=rates_hat';
Err_Q(U,:)=diag(F_hat)';
LB(U,:)=Q_bal(U,:)-abs(Err_Q(U,:));
UB(U,:)=Q_bal(U,:)+abs(Err_Q(U,:));
end
    LB_table{1,ver}={LB};
    UB_table{1,ver}={UB};

    
    
figure(ver+10)
for p=1:length(Cexp_list)
subplot(5,6,p)
plot(Q_rates(:,p),'o-b')
hold on
%errorbar(Q_bal(:,p),Err_Q(:,p),'*-r')
plot(Q_bal(:,p),'*-r')
title(Cexpl(Cexp_list(p)))

end

end
end
end
function experiment_outcomes(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2,model,EMB,MetNum,LB,UB,compSymbolList,compNameList,MB_splines,boundary_integration,Headers)

    if MB_splines == 1
     [to,P,V2,C2,SSR,MinCI2,MaxCI2,metabM,namesR,Measurements]=BDMFA_program(HPLC,Cexp,Cexpl,time,spreadsheet,page,Vars,t1,t2);
      Output.Measurements=Measurements; Output.to=to; Output.P=P; Output.V2=V2; Output.C2=C2; Output.SSR=SSR; Output.MinCI2=MinCI2; Output.MaxCI2=MaxCI2; Output.metabM=metabM; Output.namesR=namesR;
      time_exp=(time.*10)-(time(1)*10);
      h=figure(2);
      for J=1:length(metabM)
          subplot(6,7,J)
         plot(C2(:,J))
         xlim([0 600])
         hold on
         plot(time_exp(:),Measurements(:,J),'-*r')
         title(metabM(J))
      end
      figtitle=strcat(matlab.lang.makeValidName(Headers{I,1}),'_2.fig');
      savefig(h,figtitle);
%       close(h)
    end
    
    if boundary_integration == 1 
    ExometsIntegration(model,EMB,MetNum,LB,UB,compSymbolList,compNameList);   

    end        

end

