clear
close all
[~,sheet]=xlsread('C:\Users\Andre\Documents\outliers.xlsx','Sheet12');

res=zeros(38,15);
[~,Cexpl]=xlsread('C:\Users\Andre\Documents\outliers.xlsx','Sheet13');

for j = 1:11

[num,~]=xlsread('C:\Users\Andre\Documents\outliers.xlsx',sheet{j});

    for i=1:size(num,1)
    
    row=num(i,:);
    row=rmmissing(row);
    
    if row(1) == 0
        
    else
        for k=1:length(row)
        
        res(i,row(k))=res(i,row(k))+1;
        end
    end
    end
end

res=(res./11);

   figure(1)
   for x=1:38
   subplot(6,7,x)
   bar(res(x,:))
   ylim([0,1])
   title(Cexpl{x})
   end
   
res=res.^3;



  
 % agreement check
[num,~]=xlsread('C:\Users\Andre\Documents\outliers2.xlsx','Sheet1'); 

OLR_1=zeros(38,15);
OLR_1(:,:)=-1;
 for i=1:size(num,1)
    
    row=num(i,:);
    row=rmmissing(row);
    
    if row(1) == 0
        
    else
        for k=1:length(row)
        
        OLR_1(i,row(k))=1;
        end
    end
 end 
    
 val_1=res.*OLR_1;
 val_1_sum=(sum(val_1,2))/15;
 val_1_total=sum(val_1_sum);
 val_1_mean=mean(val_1_sum);
 val_1_sum(end+1)=(val_1_total);
 val_1_sum(end+1)=(val_1_mean);
 Cexpl_mean=Cexpl;
 Cexpl_mean{end+1}='Score Sum';
 Cexpl_mean{end+1}='Mean Score';
X=categorical(Cexpl_mean);
X=reordercats(X,Cexpl_mean); 
 figure(2)
   for x=1:38
   subplot(6,7,x)
   bar(val_1(x,:))
   ylim([-1,1])
   title(Cexpl{x})
   end 
   
   figure(3)
  bar(X,val_1_sum.*10)
  ylabel('Agreement Score')
  ylim([-1,1]);
  

  [num,~]=xlsread('C:\Users\Andre\Documents\outliers2.xlsx','Sheet2');
  
  OLR_2=zeros(38,15);
OLR_2(:,:)=-1;
 for i=1:size(num,1)
    
    row=num(i,:);
    row=rmmissing(row);
    
    if row(1) == 0
        
    else
        for k=1:length(row)
        
        OLR_2(i,row(k))=1;
        end
    end
 end 
 val_2=res.*OLR_2;
 val_2_sum=(sum(val_2,2))/15;
 val_2_mean=mean(val_2_sum);
 val_2_total=sum(val_2_sum);
 val_2_sum(end+1)=val_2_total;
 val_2_sum(end+1)=(val_2_mean);


 figure(4)
   for x=1:38
   subplot(6,7,x)
   bar(val_2(x,:))
%    ylim([-1,1])
   title(Cexpl{x})
   end 
   
   figure(5)
  bar(X,val_2_sum)
  ylabel('Agreement Score')
  
 
 