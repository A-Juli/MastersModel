function [Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pt_low,pt_high)

% Andrei Ligema
% 2019 DTU




% [Cexp , Cexpl, time]=exomets_data_import;

GT=0.25;

% pt_low=25; pt_high=75;
% 
for j=1:size(Cexp,2)
    interp_met=Cexp(:,j);
    missing=find(isnan(interp_met));
    if isempty(missing)==0 
    t_missing=time(missing);
    interp_met(:,2)=time;
    interp_met=rmmissing(interp_met);
    interp_fn=griddedInterpolant(interp_met(:,2),interp_met(:,1));
    interp_vals=interp_fn(t_missing);
    Cexp(missing,j)=interp_vals;
    end
     clear interp_met missing t_missing interp_fn interp_vals
end     

% figure(2)
% for j=1:size(Cexp,2)
%     clear EM RM_time ExMet MetTime
%     ExMet=Cexp(:,j);
%     MetTime=time;
%     Ex_dir=gradient(ExMet,time);
%     Ex_dir2=gradient(Ex_dir,time);
% 
%     subplot(6,7,j)
%     plot(time,Ex_dir,'-o');
%     title(Cexpl(j))
% 
% end
% 
% figure(3)
% for j=1:size(Cexp,2)
%     clear EM RM_time ExMet MetTime
%     ExMet=Cexp(:,j);
%     MetTime=time;
%     Ex_dir=gradient(ExMet,time);
%     Ex_dir2=gradient(Ex_dir,time);
% 
%     subplot(6,7,j)
%     plot(time,Ex_dir2,'-o');
%     title(Cexpl(j))
% 
% end


    



figure(10)
Clean_Met=NaN(size(Cexp,1),size(Cexp,2));
for j=1:size(Cexp,2)
    clear EM RM_time ExMet MetTime
    ExMet=Cexp(:,j);
    MetTime=time; 
    subplot(6,7,j)
    plot(time,ExMet,'-o b'),
    hold on
    ExMet(:,2)=MetTime;
    ExMet(:,3)=[1:1:size(ExMet,1)];
    
    
    flagA=0;
  
    % pass 1 for directional changes
   while flagA==0 % runs until the scan is able to reach the end of the data set without triggering a removal
       i=2;flagB=0;
             while flagB==0 % scans the current metabolite for changes in gradient indicating an outlier
    clear G0 G1 G2 G3
       G1=gradient(ExMet(i-1:i,1),ExMet(i-1:i,2));
       G2=gradient(ExMet(i:i+1,1),ExMet(i:i+1,2));
       
       if i>2 && i<(size(ExMet,1)-1)
           G0=gradient(ExMet(i-2:i-1,1),ExMet(i-2:i-1,2));
           G3=gradient(ExMet(i+1:i+2,1),ExMet(i+1:i+2,2));
           
       end
             
              
       if (G1(1)*G2(1))<0

           if i==2 || i==(size(ExMet,1)-1)
           ExMet(i,1)=NaN;
           ExMet=rmmissing(ExMet);
           flagB=1;    
           else
           
           if (G0(1)*G1(1))>=0 && (G2(1)*G3(1))>=0
           else
           ExMet(i,1)=NaN;
           ExMet=rmmissing(ExMet);
           flagB=1;       
           end
           end
       end
          
       i=i+1;  

       if i==(size(ExMet,1))
           flagA=1;
           flagB=1;
       end  
            end
   
%    for k=1:size(ExMet,1)
%        Clean_Met(ExMet(k,3),j)=ExMet(k,1);
%    end
%    
   
   end
   
   ExMet(1:size(ExMet,1),4)= zeros;
   
   J=(size(ExMet,1)-2);
   % pass 2 for smoothing of curves
   for I=2:J
       clear G0 G1 G2 G3
   G1=gradient(ExMet(I:I+1,1),ExMet(I:I+1,2));
   G2=gradient(ExMet(I+1:I+2,1),ExMet(I+1:I+2,2));
   G3=gradient(ExMet(I:I+2,1),ExMet(I:I+2,2));
   if abs(G3(1))<abs(G1(1)*GT) || abs(G3(1))<abs(G2(1)*GT) 
           if ExMet(I,4)~=1
       ExMet(I,1)=NaN;
       if I~=1
       ExMet(I+1,4)=1;   ExMet(I-1,4)=1; 
       else
       ExMet(I+1,4)=1;
       end
           end
   end
   
   end
   for k=1:size(ExMet,1)
       Clean_Met(ExMet(k,3),j)=ExMet(k,1);
   end

   ExMet=rmmissing(ExMet);
   plot(ExMet(:,2),ExMet(:,1),'-* r')
   title(Cexpl(j))
%    hold on
%    plot(time,TF,'*')
%   

    end

end

