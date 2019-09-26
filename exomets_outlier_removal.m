function [Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pt_low,pt_high)

% Andrei Ligema
% 2019 DTU




% [Cexp , Cexpl, time]=exomets_data_import;


% pt_low=25; pt_high=75;



figure(2)
for j=1:size(Cexp,2)
    clear EM RM_time ExMet MetTime
    ExMet=Cexp(:,j);
    MetTime=time;
    Ex_dir=gradient(ExMet,time);

%     subplot(6,6,j)
%     plot(time,Ex_dir,'-o');
%     title(Cexpl(j))

end

figure(10)

for j=1:size(Cexp,2)
    clear EM RM_time ExMet MetTime
    ExMet=Cexp(:,j);
    MetTime=time; 
    subplot(6,6,j)
    plot(time,ExMet,'-o b'),
    hold on
    Ex_dir=gradient(ExMet,time);
    [~, TF]=rmoutliers(Ex_dir,'percentile',[pt_low pt_high]);
   
    
    
   for i=1:length(TF)
       if TF(i)==1
%            MetTime(i)=NaN;
%            MetTime(i+1)=NaN;
           ExMet(i)=NaN;
%            ExMet(i+1)=NaN;
           
       else
       end
      Clean_Met(i,j)=ExMet(i); 
   end
   ExMet(:,2)=MetTime;
   
   ExMet=rmmissing(ExMet);
   plot(ExMet(:,2),ExMet(:,1),'-* r')
   title(Cexpl(j))
%    hold on
%    plot(time,TF,'*')
%   

end

end

