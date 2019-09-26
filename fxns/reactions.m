%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% to generate stoichiomatric matrix from a list of reactions
function [M,Mc,compr,compp,names,namesc,transport]=reactions(Reactions)
reaction=Reactions(:,2);
reaction(1)=[];
rea=cell(size(reaction,1),36);
re=cell(size(reaction,1),36);
srea=zeros(size(reaction,1),36);
prd=cell(size(reaction,1),36);
pr=cell(size(reaction,1),36);
sprd=zeros(size(reaction,1),36);
compr=cell(size(reaction,1),36);
compp=cell(size(reaction,1),36);
transport=ones(size(reaction,1),1);
for i= 1:size(reaction,1)
    rxn=reaction{i};
    c=findstr('[',rxn);
    if c(1)==1
        transport(i)=0;
        compr{i,1}=rxn(2:2);
        compp{i,1}=rxn(2:2);
        rxn=rxn(4:end);
        [reactants,products]=divideReaPrd(rxn);
        [re(i,:),srea(i,:)]=divideMetab(reactants);
        [pr(i,:),sprd(i,:)]=divideMetab(products);
        rea(i,:)=AddComp(re(i,:),compr{i,1});
        prd(i,:)=AddComp(pr(i,:),compp{i,1});
    else
        [reactants,products]=divideReaPrd(rxn);
        [rea(i,:),srea(i,:)]=divideMetab(reactants);
        [prd(i,:),sprd(i,:)]=divideMetab(products);
        [re(i,:),compr(i,:)]=RemoveComp(rea(i,:));
        [pr(i,:),compp(i,:)]=RemoveComp(prd(i,:));
    end  
  %fprintf('\n %1.0f.\n',i);
end
[names,M]=CreateM(re,pr,srea,sprd);
[namesc,Mc]=CreateM(rea,prd,srea,sprd);
    

