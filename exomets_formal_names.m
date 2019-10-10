function [EMB,MetNum,model]=exomets_formal_names (model,Cexpl)

modelNew=model;
EMB=Cexpl(5:34);
compSymbolList=modelNew.comps;
compNameList=modelNew.compNames;

% search for formal metabolite names in model
EMB_rem=[];
for i=1:length(EMB)
    clear MetLocs MetIdx MetIdx_ex MetIdx_ex_loc MetCode
    MetLocs=strfind(lower(model.metNames),lower(EMB(i)));
    MetIdx=find(~cellfun('isempty',MetLocs)); % returns list of indexes of matches
    MetIdx_ex=strfind(model.mets(MetIdx),'[e]'); 
    MetIdx_ex_loc=find(~cellfun('isempty',MetIdx_ex));
    
    MetCode=model.mets(MetIdx(MetIdx_ex_loc));
    
    if isempty(MetCode)==1
        fprintf('No metabolite matches found for %s \n',(EMB{i}));
        found=model.mets(MetIdx)
        new_met=input('Are the idenfitied metabolites, do you wish to add a new metabolite? \n','s');
        if new_met=='Y'
            metID=input('Input MetID:\n','s');
            metName=input('Input MetName:\n','s');
            formula=input('Input MetFormula:\n','s');
            Charge=input('Input Charge:\n');
            modelNew=addMetabolite(modelNew,metID,'metName',metName,'metFormula',formula,'Charge',Charge);
            MetCode={metID};
            Met_no=(length(modelNew.mets))+1;
        else
            rem_met=input('Do you wish to ignore this metabolite? \n','s');
            if rem_met=='Y'
                EMB_rem(i)=i;
                MetCode=EMB(i);
            else
                error('aborting process');
            end
        end
    else
        if size(MetCode,1)>=2
        fprintf('Multiple matches found for metabolite %s \n',EMB{i});
        MetCode 
        metpick=input('\n Please input the number of the appropriate metabolite: \n');
        MetCode=MetCode(metpick);
        Met_no=MetIdx(MetIdx_ex_loc(metpick));  
        else 
        Met_no=MetIdx(MetIdx_ex_loc);    
        end    
        
    end           
    

    
    MetNum(i)=Met_no;
    EMB(i)=MetCode;
end
if isempty(EMB_rem)~=1
EMB_rem=rmmissing(EMB_rem);
EMB(:,EMB_rem)=[];    
end
model=modelNew;
end