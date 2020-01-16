function [dataOUT,rates_hat,F_hat]=GMEA(model,Q_rates,se,MetNum,EMB,U,Biom)


%gross measurement error analysis (GMEA)
%use hypothesis testing to detect gross error in measurement set
%author Lake-Ee Quek, The University of Queensland, AIBN
%02-Oct-2008

% Modified version for exomets_process QMCM DTU Biosustain
% Andrei Ligema 2019

% Inputs    [model] GSM, no boundaries or objective required
%           [Q_rates] Imported Flux across [k] clusters, developed by
%           k-means procedure. May be run based on a single replicate
%           [se] table of standard error for each measured metabolite
%           [MetNum] Index of metabolite locations linking Q_rates and
%           model.mets
%           [EMB] List of metabolite names validated by
%           exomets_formal_names
%           [U] tracks experimental replicates
%           [Biom] unused
% Outputs   [dataOUT] Compiled structure of output variables
%           [rates_hat] balanced flux predictions
%           [F_hat] Covariance matrix based on balanced predictions

% Methods from the original code have been left intact and should function
% if desired. Their purpose was to generate every possible iteration of the
% reduced redundancy matrix (n choose k, where n is the rows in rm and k is
% the rank of rm) and then test those combinations to find the version with
% the lowest condition number. While this may generate the optimal version
% of the reduced redundancy matrix it is VERY computationally intensive to
% perform with even a reduced GSM (order of 370C7), only use if generating
% the best conditioned Rred is vital. Otherwise the replacement code should
% produce a good enough Rred almost immediately.




% Depricated, part of the original code
%select reaction to delete
% rxn_del ={'EF0031'};
% for i = 1:length(rxn_del)
%     rxnID = rxn_del{i};
%     for j = 1:size(rxnTable,1)
%         rxnString = rxnTable{j,1};
%         if strcmp(rxnString,rxnID)
%             rxnTable(j,:)=[];
%             stoicMat(j,:)=[];
%             break
%         end
%     end
% end
%read input data
%column 1: rxnMea (string);%column 2: rates (double)
%column 3: standard error (double)
% disp('input file must have 3 columns: rxnID, rates, standard error')
% ratesFile = input('specify name of file containing measured rates  >','s');
% [rxnMea,rates,se] = textread(ratesFile,'%s%f%f');
%se=se;


rxnTable=model.rxns;
stoicMat=model.S;

rates=Q_rates(U,:)';
% locate measured exchange reactions 
[Ex]=rxn_loc(model,MetNum,EMB);
Ex=Ex(:,2);

clear dataOUT
rxnMea=model.rxns(Ex);
%separate stoic matrix (Sm'vm + Sc'vc = 0)
%into measured and calculated
rowMea=[];
for i = 1:length(rxnMea)
    rxnMeaID = rxnMea{i};
    for j = 1:size(rxnTable,1)
        rxnID = rxnTable{j,1};
        if strcmp(rxnMeaID,rxnID)
            rowMea(i)= j;
        end
    end
end
rowCalc=[];
for i = 1:size(stoicMat,1)
    if length(find(rowMea==i)) == 0
        rowCalc=[rowCalc i];
    end
end
stoicMea = stoicMat(rowMea,:);
stoicCalc = stoicMat;
stoicCalc(rowMea,:)=[];

%calculate the redundancy matrix
rm = stoicMea' - stoicCalc'* pinv(stoicCalc')* stoicMea';

%rank shows degree of redundancy
disp(['rank:   ' num2str(rank(rm,1e-3))]);

% pre-treating the redundancy matrix
% locs=[1:size(rm,1)];
% locating empty rows
rm_sum=sum(abs(rm),2);
z_locs=find(rm_sum<1E-5);
rm(z_locs,:)=[];
%locating duplicates
I=1;
while 1
   
ident_test=(abs(rm)-abs(rm(I,:)));
test_rows=find(abs(sum(ident_test,2))<1E-5);
if length(test_rows)>1
    rm(test_rows,:)=[];
    I=1;
else
    I=I+1;
end

if I==(size(rm,1)+1)
    break
end
    
    
end


%rank should be 2 or greater, pick n rank rows for metabolite balance
%find the combination of rows that give the smallest conditional number
%generate combination
count = rank(rm,1e-3);

% Original Code shown for comparison, this procedure would generate all
% combinations for N Choose K in an iterative fashion, then check the
% condition number for each composition by SVD. Needless to say this is
% very inefficient method with models more significant than toy tests.
% DO NOT USE, ever, a waste of computation

% maxCount = size(rm,1);
% combStore=[];
% firstComb = [1:count];
% lastComb =[maxCount-count+1:maxCount];
% combStore=[combStore;firstComb];
% while 1
%     latestComb = combStore(end,:);
%     updatedCell = 0;
%     for cellCheck = count:-1:1
%         if latestComb(cellCheck) < lastComb(cellCheck)
%             latestComb(cellCheck) = latestComb(cellCheck) + 1;
%             updatedCell = cellCheck;
%             break;
%         end
%     end
%     if updatedCell > 0
%         increment = 1;
%         for cellMod = updatedCell+1:count
%             latestComb(cellMod) = latestComb(updatedCell)+increment;
%             increment = increment + 1;
%         end
%         combStore=[combStore;latestComb];
%     end
%     if cellCheck == 1 && updatedCell == 0
%         break
%     end
% end
% 

% A much faster way of generating the iterations of N Choose K! Depricated
% with the methods used to generate the reduced redundancy matrix
% introduced below
% v=1:1:(size(rm,1));
% combStore=nchoosek(v,count);
% 

% New method
% Pre-treating the Redundancy Matrix, removing rows with 0 information and
% rows that show obvious linear dependence due to identical coefficients
% Remaining rows are then sorted based on how close the (abs) values of their
% non-zero elements are to unity. As these rows are likely to contribute to
% a better conditioned Reduced Redundancy matrix
for i=1:size(rm,1)   
    loc=find(abs(rm(i,:))>1e-8);
    U(i,1)=i;    
    U(i,2)=sum(abs(rm(i,loc)));   
    U(i,3)=length(loc);
    U(i,4)=U(i,2)/U(i,3);
end

U=sortrows(U,4,'descend');

% generates a reduced redundancy matrix based on the above criteria 
% starting with the top it attempts to find a set with full rank. 
% This procedure needs to be updated in order to be able to discard the
% first row and repeat the procedure if no solution is found including row
% 1

T=2;
P=2;
rowSelect=[1];
% T tracks the number of entries in rowselect
% P tracks the current inclusion candidate
while 1
 % generate combination
    rowSelect(T)=P;
    % Test combination
    
    if rank(rm(rowSelect,:),1e-3)==T
        % Passed Rank test increment T and P
        T=T+1;
        P=P+1;
        if T==count+1
            break
        end
    else
        % Rank test failed, try the next position
        P=P+1;
    end
     
    
end



% 
% Depricated, associated with original code. Use in combination with
% previously commented code to create a fully optimal Rred. But bear in
% mind that even with a pre-treated redundancy matrix this will likely
% require millions of Single Value Decompositions. There is probably
% something better you can be doing with your time.

% %screen combination
% rowSelect =combStore(1,:);
% condStore=cond(rm(rowSelect,:));
% for i = 1:size(combStore,1)
%     if cond(rm(combStore(i,:),:)) < condStore
%         condStore = cond(rm(combStore(i,:),:));
%         rowSelect=combStore(i,:);
%     end
% end

%generate reduced redundancy matrix
rm_reduced = rm(rowSelect,:);
%calculate residual error
res_error = rm_reduced*rates;
%generate variance-covariance matrix
%assume independent measurements, error occupies diagonal elements
F=diag(se.^2); 
P=rm_reduced*F*rm_reduced';

%test hypothesis, DoF 2, % s chi value of    5.99
h_test = res_error'*inv(P)*res_error;
fprintf('test score:\t%1.3f\n',h_test);
fprintf('chi score cut-off:\t%1.3f\n',chi2inv(0.95,count));

%minimize error and adjust measured rates
del_hat = F*rm_reduced'*inv(P)*rm_reduced*rates;
rates_hat = rates - del_hat;

%calculate variance of adjusted rates
F_hat = F - F*rm_reduced'*inv(P)*rm_reduced*F;

%calculate unknown rates and associated variance
rates_calc_hat = - pinv(stoicCalc')*stoicMea'*rates_hat;
F_calc_hat = pinv(stoicCalc')*stoicMea'*F_hat*stoicMea*pinv(stoicCalc')';

dataOUT = struct;
dataIN = [rates_hat F_hat];
dataOUT.rates_hat=rates_hat; dataOUT.del_hat=del_hat; dataOUT.F_hat=F_hat;
dataOUT.rates_calc_hat=rates_calc_hat; dataOUT.F_calc_hat=F_calc_hat;

end

function [Ex]=rxn_loc(model,MetNum,EMB)

% Search for corresponding exchange reactions in the model, a reduced clone
% of exomets_integration procedure for identifying exchange reactions based on
% metabolites
rxns=cell(length(EMB),1);
Ex=zeros(length(EMB),2);

% [EMB,MetNum,model]=exomets_formal_names (model,Cexpl);



for i=1:(length(MetNum))
        RxL=find(abs(model.S(MetNum(i),:))==1);
        rxns{i}=RxL;
    for j=1:(length(rxns{i}))
        clear RL MR
        RL=rxns{i}(1,j);
        MR=find(model.S(:,RL));
        
             RxnNet(i,j)=strcmp(model.metNames{MR(1),1},model.metNames{MR(2),1});
     
  %        RxnNet(i,j)=sum(abs(model.S(:,rxns{i}(1,j))));
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

ExT(:,3)=model.rxnNames(ExT{:,2});

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
ExT(:,3)=array2table(model.rxns(Ex(:,2)));


end