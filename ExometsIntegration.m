%Import data from spreadsheets, used for testing/
ExoMetsList=readtable('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\Exomets.xlsx','Sheet','MetsList');
ExoMetsList=table2struct(ExoMetsList,'ToScalar',true);
ExoMetsTabP71=readtable('C:\Users\anjuli\Documents\MATLAB\cobratoolbox\CHO\Exomets.xlsx','Sheet','P71');
ExoMetsP71=table2struct(ExoMetsTabP71,'ToScalar',true);
LB(1:length(ExoMetsList.Mets))=1e-7;
UB(1:length(ExoMetsList.Mets))=900;
ExoMetsList.LB=LB';
ExoMetsList.UB=UB';
%%%%%%%%%%%%%%% ABOVE SEGMENT IMPORTED FOR VALIDATION

% Input format: 
% one table to be converted to a structure containing
% Col 1: metabolite codes ie. glc_D_[e]
% Col 2: corresponding row entry in the mets list for the model
% Col 3: Lower boundaries for the target exchange reactions
% Col 4: Upper boundaries for the target exchange reactions



% Exomets data integration
% Inputs: matrix containing metabolite numbers for relevant exomets and
% their derived upper and lower boundaries
% & the COBRA model for the target organism.
% Outputs: Updated COBRA model incorporating the exometabolomic data as
% constraints for exchange reactions.
clear rxns RxnNet Ex ExR RxL FNum FNam FNamT FNumT Fails
model=model_iCHOv1_S_Base;
modelNew=model;
EMB=ExoMetsList;
compSymbolList=model.comps;
compNameList=model.compNames;

% Search for corresponding exchange reactions in the model
rxns=cell(length(EMB.Mets),1);
Ex=zeros(length(EMB.Mets),2);


for i=1:(length(EMB.MetNo))
        RxL=find(abs(model.S(EMB.MetNo(i),:))==1);
        rxns{i}=RxL;
    for j=1:(length(rxns{i}))
        RxnNet(i,j)=sum(abs(model.S(:,rxns{i}(1,j))));
    end
    %searches for an exchange reaction based on the absolute sum of
    %stoichometries which should be 1 for an exchange reaction only
    ExR=(find(RxnNet(i,:)==1));
    if isempty(ExR)==0
        Ex(i,1)=EMB.MetNo(i);
        Ex(i,2)=rxns{i}(1,ExR);
        
    else
        % if no reaction satisfies the conditions for Exchange the entry is
        % given a flag value to prompt further investigation
        ExR=NaN;
        Ex(i,1)=EMB.MetNo(i);
        Ex(i,2)=ExR;
        
    end
end
clear i j
% Produces the matrix Ex which contains the rows and columns of the
% exchange reactions for the target metabolites in the S matrix.
% Converting to a table to import the reaction names
ExTNam=["Met_Index","Ex_Rxn","Rxn_Name","LB","UB"];
ExTTyp=["double","double","string","double","double"];
ExT=table('size',[length(EMB.Mets),5],'VariableNames',ExTNam,...
    'VariableTypes',ExTTyp);
ExTP=array2table(Ex);
ExT(:,1:2)=ExTP(:,1:2);
ExBC(:,1)=array2table(EMB.LB);
ExBC(:,2)=array2table(EMB.UB);
ExT(:,4:5)=ExBC(:,1:2);
clear ExTP ExTNam ExTTyp ExBC

% export list of failures
if isempty(find(isnan(Ex)))==0
      
FNum=EMB.MetNo((find(isnan(Ex))-length(Ex)));
FNumT=array2table(FNum);
FNam=EMB.Mets((find(isnan(Ex))-length(Ex)));
FNamT=cell2table(FNam);
FLBs=EMB.LB((find(isnan(Ex))-length(Ex)));
FLBs=array2table(FLBs);
FUBs=EMB.UB((find(isnan(Ex))-length(Ex)));
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
ExT(:,3)=array2table(model.rxns(Ex(:,2)));

% Perform boundary integration with the COBRA model   
% For successfully isolated reactions
   for i=1:height(ExT)
       RxnNam=char(ExT{i,3});
       modelNew=changeRxnBounds(modelNew,RxnNam,ExT{i,4},'l');
       modelNew=changeRxnBounds(modelNew,RxnNam,ExT{i,5},'u');
   end
   % Create new reactions from those without exchanges
   for j=1:height(Fails)
       MetID=char(model.mets(Fails{j,1}));
       NewRxnID=length(modelNew.rxns)+1;
       modelNew=addExchangeRxn(modelNew,MetID,Fails{j,3},Fails{j,4});
       RID='EXCHANGE';
       modelNew.subSystems{NewRxnID}={RID};
              
   end
   modelNew.subSystems=convertStringsToChars(modelNew.subSystems);
   clear RID

   % Validate by printing constraints
   clear i j
   printConstraints(modelNew,-1000,1000)
   fprintf('new model will be saved in SBML .xml format');
SaveState=input('Do you want to save the new model Y/N? ','s');
if SaveState=='Y'
    ModelNam=input('please enter the filename for the updated model: \n','s');
    ModelNam=char(ModelNam);
sbmlModel=writeSBML(modelNew,ModelNam,compSymbolList,compNameList);
else
end
fprintf('cleaning up \n');

clear sbmlModel
clear model
clear ExoMetsList ExoMetsTabP71 ExoMetsP71 LB UB compSymbolList compNameList
clear EMB rxns Ex RxL rxns RxnNet i j ExR ExT Fails RxnNam MetID SaveState
clear ModelNam NewRxnID
    