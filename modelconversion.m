    

modelDirectory='C:\Users\anjuli\Documents\MATLAB\cobratoolbox\ExoMets';
           modelFileName_base='Ivan_mm.xml';
           modelFileName_base = [modelDirectory filesep modelFileName_base];
           model = readCbModel(modelFileName_base);
           compSymbolList=model.comps;
           compNameList=model.compNames;
       modelDirectory='C:\Users\anjuli\Documents\MATLAB\cobratoolbox\ExoMets';
       modelFileName_base='GSM_mouse_trimmed_BDMFA.xml';
       modelFileName_base = [modelDirectory filesep modelFileName_base];
       model_base = readCbModel(modelFileName_base);


model=removeRxns(model,{'BIOM_AA','BIOM_DNA','BIOM_RNA','BIOM_LIP','BIOM_CARBO','BIOM_T'});

rxn_label='BIOToT';

Mets=model_base.metNames(find(model_base.S(:,262)));
Stoich=model_base.S(find(model_base.S(:,262)),262);

for i=1:length(Mets)
    
Metdirty=Mets{i};
Metdirty=strrep(Metdirty,'-','');
Metdirty=strrep(Metdirty,'(','');
Metdirty=strrep(Metdirty,')','');
match_mets=strmatch(Metdirty,model.metNames);

% if isempty(match_mets)==1
%    fprintf('%s \n',Mets{i})
%    metname=input('No match found, please input the name manually for \n','s');
%    match_mets=strmatch(metname,model.metNames);    
% end

if match_mets >= 2
   loc=contains(model.mets(match_mets),'cytosol');       
   match_mets=match_mets(loc); 
end

list(i)=match_mets;


end
list=list';
list=model.mets(list);
model=addReaction(model,rxn_label,'metaboliteList',list,'stoichCoeffList',Stoich,'reversible',false);
 model.c(1:length(model.rxns),1)=0;
   model=changeObjective(model,'BIOToT',1);
sbmlModel=writeSBML(model,'Ivan_mm_biom.xml',compSymbolList,compNameList);