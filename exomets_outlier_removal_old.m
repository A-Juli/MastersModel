function [Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pt_low,pt_high)

% Andrei Ligema
% 2019 DTU




% [Cexp , Cexpl, time]=exomets_data_import;


% pt_low=25; pt_high=75;

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

figure(2)
for j=1:size(Cexp,2)
    clear EM RM_time ExMet MetTime
    ExMet=Cexp(:,j);
    MetTime=time;
    Ex_dir=gradient(ExMet,time);
    Ex_dir2=gradient(Ex_dir,time);

    subplot(6,7,j)
    plot(time,Ex_dir,'-o');
    title(Cexpl(j))

end

figure(3)
for j=1:size(Cexp,2)
    clear EM RM_time ExMet MetTime
    ExMet=Cexp(:,j);
    MetTime=time;
    Ex_dir=gradient(ExMet,time);
    Ex_dir2=gradient(Ex_dir,time);

    subplot(6,7,j)
    plot(time,Ex_dir2,'-o');
    title(Cexpl(j))

end


    



figure(10)

for j=1:size(Cexp,2)
    clear EM RM_time ExMet MetTime
    ExMet=Cexp(:,j);
    MetTime=time; 
    subplot(6,7,j)
    plot(time,ExMet,'-o b'),
    hold on
    Ex_dir=gradient(ExMet,time);
    Ex_dir=gradient(Ex_dir,time);
    [~, TF]=rmoutliers(Ex_dir,'percentile',[pt_low pt_high]);
   
    
    Clean_Met(1,j)=ExMet(1);  
   for i=2:(length(TF)-1)
       if TF(i)==1
       G1=gradient(Cexp(i-1:i,j),time(i-1:i));
       G2=gradient(Cexp(i:i+1,j),time(i:i+1));
       if (G1(1)*G2(1))<0
           ExMet(i)=NaN;
       end
       
       else
       end
      Clean_Met(i,j)=ExMet(i); 
   end
   Clean_Met(i+1,j)=ExMet(i+1);  
   ExMet(:,2)=MetTime;
   
   ExMet=rmmissing(ExMet);
   plot(ExMet(:,2),ExMet(:,1),'-* r')
   title(Cexpl(j))
%    hold on
%    plot(time,TF,'*')
%   

end

end

