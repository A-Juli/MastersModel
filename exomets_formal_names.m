function [EMB,MetNum,model,Cexp_list]=exomets_formal_names (model,Cexpl,ExTag,Util_mets)

% Module designed to handle inconsistencies between labelling of
% metabolites in Cexpl (assuming this was imported from the results sheet)
% and their naming convention in the selected model. Returns a list of
% numerical indices for the measured metabolites, simplifying integration
% and interrogation of the model
% If a metabolite is found to be missing this function is able to add it to
% the model and produce an updated version to be used by future replicates
% of this experiment 
% Inputs:   [model] A Constraint Based GSM
%           [Cexpl] Cell list of metabolites measured
%           [ExTag] The characters used in your GSM to identify metabolites
%           in the external compartment in model.mets
%           [Util_mets] Index of metabolites to be considered, used if more
%           metabolites have been measured than are implemented in the
%           model
% Note: due to inconsistencies in metabolite labelling between GSMs the
% user may need to adjust the search or display fields to make proper use
% of this function. Relevant lines will be highlighted.
% 
% By sharing index values between EMB, MetNum and Cexp_list we tie together
% the names and locations in data and model for utilised metabolites,
% making it less likely to get tangled in the future!
% Outputs:  [EMB] Cell list of metabolites which are participating in the
%           model, drawn from model.mets to ensure character matches
%           [MetNum] List of linear indices within model.mets of utilised
%           metabolites
%           [Cexp_list] list of indexes within Cexp that associate with
%           utilised metabolites
%           [model] updated version of the GSM including newly written
%           reactions

modelNew=model; %clones the import model to allow for modifications

if isempty(Util_mets)==1
    warning('List of used metabolites is empty, presuming at all measured metabolites are intended for inclusion')
    EMB=Cexpl;
    Cexp_list=[1:length(EMB)];
else    
    EMB=Cexpl(Util_mets); % imports the "raw" list from input
    Cexp_list=Util_mets;
end

% compSymbolList=modelNew.comps;
% compNameList=modelNew.compNames;

ignoremissing=0; 
% debug variable used to automatically ignore missing metabolites which
% would otherwise require a new addition or throw an error. Most certainly
% not reccomended for live use

% search for formal metabolite names in model
EMB_rem=[]; % holds list of failures to remove, only becomes active if a missing metabolite is ignored
for i=1:length(EMB)
    clear MetLocs MetIdx MetIdx_ex MetIdx_ex_loc MetCode
    MetLocs=strfind(lower(model.mets),lower(EMB(i))); % locates matches for metabolite, this searches metNames which is usually given in natural language and so easier to search
    MetIdx=find(~cellfun('isempty',MetLocs)); % returns list of indexes of matches
    MetIdx_ex=strfind(model.mets(MetIdx),ExTag); % checks the candidate list MetIdx in model.mets to identify metabolites flagged as being in the external compartment
    MetIdx_ex_loc=find(~cellfun('isempty',MetIdx_ex)); 
    
    MetCode=model.mets(MetIdx(MetIdx_ex_loc));
    
    if isempty(MetCode)==1 % If no metabolites returned solid matches in both identifications you will see this
        if ignoremissing==0        
        fprintf('No metabolite matches found for %s \n',(EMB{i})); 
        found=model.mets(MetIdx) 
        warning('If this list contains clearly correct metabolites then you should check variable [ExTag] or exomets_formal_names search fields ~line 50')
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
            rem_met=input('Do you wish to ignore this metabolite? (not reccomended) \n','s');
            if rem_met=='Y'
                EMB_rem(i)=i;
                MetCode=EMB(i);
            else
                error('aborting process');
            end
        end
        else
              EMB_rem(i)=i;
              MetCode=EMB(i);    
        end % end of match failure method
    else 
        if size(MetCode,1)>=2 % multiple potential matches found, this may be because of shared strings such as Isoleucine / leucine, user should select the appropriate metabolite
            % if this becomes annoying consider implementing stricter
            % matching criteria, not reccomended for regular use as it
            % depends on the model not using crazy naming schemes
        fprintf('Multiple matches found for metabolite %s \n',EMB{i}); 
        found=model.mets(MetIdx(MetIdx_ex_loc))
        metpick=input('\n Please input the number of the appropriate metabolite: \n');
        if isnumeric(metpick)~=1
            error('Invalid specification, Aborting')
        end
        MetCode=model.metNames(MetIdx(MetIdx_ex_loc(metpick)));
        Met_no=MetIdx(MetIdx_ex_loc(metpick));  
        else % one match found
        Met_no=MetIdx(MetIdx_ex_loc);    
        end    
        
    end
    
    
    
    MetNum(i)=Met_no; % for export
    EMB(i)=MetCode; % replaces the natural language metabolite name with the formal code from model.mets
end

% If metabolites were ignored they are purged from the list here
if isempty(EMB_rem)~=1
EMB_rem=rmmissing(EMB_rem);
EMB_rem=EMB_rem(EMB_rem~=0);
EMB(:,EMB_rem)=[];   
MetNum(:,EMB_rem)=[];
Cexp_list(:,EMB_rem)=[];
end

% if the model was updated then the new version will be passed to the main
% procedure
model=modelNew;
end