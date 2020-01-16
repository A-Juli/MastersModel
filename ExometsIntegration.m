  function ExometsIntegration(model,EMB,MetNum,LB,UB,compSymbolList,compNameList)
% Andrei Ligema DTU/CFB 2019

% Exomets data integration
% Inputs: matrix containing metabolite numbers for relevant exomets and
% their derived upper and lower boundaries
% & the COBRA model for the target organism.
% Outputs: Updated COBRA model incorporating the exometabolomic data as
% constraints for exchange reactions.

% Procedure: using the formalised metabolite positions from
% exomets_formal_names, identifies the specific exchange reactions involved
% in transporting the target exometabolites and sets the boundary
% constraints based on an input set (this may be from any procedure). 
% It then sets an objective function and outputs an updated model in SBML
% format for further use. WARNING: There appear to be some issues with the
% function sbmlwrite, the models produced may not function correctly. 
clear rxns RxnNet Ex ExR RxL FNum FNam FNamT FNumT Fails





ExoMetsList=EMB;


for k=1:size(LB,1)
lb=LB(k,:);
ub=UB(k,:);



% Search for corresponding exchange reactions in the model
rxns=cell(length(EMB),1);
Ex=zeros(length(EMB),2);



modelNew=model;

% open all drains
modelNew.lb(strmatch('EX_',modelNew.rxns))=-1000;
modelNew.ub(strmatch('EX_',modelNew.rxns))=1000;

for i=1:(length(MetNum))
        RxL=find(abs(modelNew.S(MetNum(i),:))==1);
        rxns{i}=RxL;
    for j=1:(length(rxns{i}))
        clear RL MR
        RL=rxns{i}(1,j);
%         MR=find(modelNew.S(:,RL));
        
       
%             RxnNet(i,j)=strcmp(modelNew.metNames{MR(1),1},modelNew.metNames{MR(2),1});
     
         RxnNet(i,j)=sum(abs(modelNew.S(:,rxns{i}(1,j))));
    end
    %searches for an exchange reaction based on the absolute sum of
    %stoichometries which should be 1 for an exchange reaction only
    ExR=(find(RxnNet(i,:)==1));
    
    if isempty(ExR)==0
        Ex(i,1)=MetNum(i);
        Ex(i,2)=rxns{i}(1,ExR);
        
    else
        % if no reaction satisfies the conditions for Exchange the entry is
        % given a flag value to prompt further investigation
        ExR=NaN;
        Ex(i,1)=MetNum(i);
        Ex(i,2)=ExR;
        
    end
end
clear i j
% Produces the matrix Ex which contains the rows and columns of the
% exchange reactions for the target metabolites in the S matrix.
% Converting to a table to import the reaction names
ExTNam=["Met_Index","Ex_Rxn","Rxn_Name","LB","UB"];
ExTTyp=["double","double","string","double","double"];
ExT=table('size',[length(EMB),5],'VariableNames',ExTNam,'VariableTypes',ExTTyp);
ExTP=array2table(Ex);
ExT(:,1:2)=ExTP(:,1:2);
ExBC(:,1)=array2table(lb');
ExBC(:,2)=array2table(ub');
ExT(:,3)=model.rxnNames(ExT{:,2});
ExT(:,4:5)=ExBC(:,1:2);
clear ExTP ExTNam ExTTyp ExBC

% export list of failures
if isempty(find(isnan(Ex)))==0
      
FNum=MetNum((find(isnan(Ex))-length(Ex)));
FNumT=array2table(FNum);
FNam=EMB((find(isnan(Ex))-length(Ex)));
FNamT=cell2table(FNam);
FLBs=lb((find(isnan(Ex))-length(Ex)));
FLBs=array2table(FLBs);
FUBs=ub((find(isnan(Ex))-length(Ex)));
FUBs=array2table(FUBs);
FailTabNames=["Met_Index","Met_Name","LB","UB"];
FailTabTypes=["double","string","double","double"];
Fails=table('size',[length(FNum),4],'VariableNames',FailTabNames,...
    'VariableTypes',FailTabTypes);

Fails(:,1)=FNumT;
Fails(:,2)=FNamT;
Fails(:,3)=FLBs;
Fails(:,4)=FUBs;
Fails
fprintf('The above metabolites returned no exchange reactions based on stoichometry \n');
fprintf('New Reactions will be added based on the name of the metabolite \n');
clear FNum FNumT FNam FNamT FLBs FUBs
clear FailTabNames FailTabTypes
else
    fprintf('All metabolites located exchange reactions based on stoichometry \n');
end

% clean failures from ExT and locate reaction names for COBRA function
ExT=rmmissing(ExT,'MinNumMissing',2);
Ex=rmmissing(Ex);
ExT(:,3)=array2table(modelNew.rxns(Ex(:,2)));

% [result]=boundary_sequential_check(ExT,model,rxns);

% Perform boundary integration with the COBRA model   
% For successfully isolated reactions
   for i=1:height(ExT)
       RxnNam=char(ExT{i,3});
       modelNew=changeRxnBounds(modelNew,RxnNam,ExT{i,4},'l');
       modelNew=changeRxnBounds(modelNew,RxnNam,ExT{i,5},'u');
   end
   % Create new reactions from those without exchanges
 if isempty(find(isnan(Ex)))==0  
   for j=1:height(Fails)
       MetID=char(modelNew.mets(Fails{j,1}));
       NewRxnID=length(modelNew.rxns)+1;
       modelNew=addExchangeRxn(modelNew,MetID,Fails{j,3},Fails{j,4});
       RID='EXCHANGE';
       modelNew.subSystems{NewRxnID}={RID};
              
   end
   
   modelNew.subSystems=convertStringsToChars(modelNew.subSystems);
 end
   clear RID

   % add objective function
   modelNew.c(1:length(rxns),1)=0;
   modelNew=changeObjective(modelNew,'biomass_cho_producing',1);
   
   % Validate by printing constraints
   clear i j
   printConstraints(modelNew,-1000,1000)
   fprintf('new model will be saved in SBML .xml format');
SaveState=input('Do you want to save the new model Y/N? ','s');
if SaveState=='Y'
    ModelNam=input('please enter the filename for the updated model: \n','s');
    ModelNam=char(ModelNam);
    ModelList{k}=ModelNam;
sbmlModel=writeSBML(modelNew,ModelNam,compSymbolList,compNameList);
else
end
fprintf('cleaning up \n');

clear sbmlModel

clear ExoMetsList ExoMetsTabP71 ExoMetsP71 lb ub 
clear  rxns Ex RxL rxns RxnNet i j ExR ExT Fails RxnNam MetID SaveState
clear ModelNam NewRxnID
    
end
  end