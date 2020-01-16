function [Clean_Met]=exomets_outlier_removal(Cexp,Cexpl,time,pre_gradient,pre_gradient_range,exp_weight,Exp_Err)

% Andrei Ligema
% 2019 DTU/CFB

% A procedure to filter short-length metabolic data with no or few repeats.
% Functions by using an exponentially weighted moving mean of gradient to
% test the predictability of subsequent datapoints based on the previous
% curvature within experimental error. Points that fall out of this range
% are flagged for removal. 
% Flagged points are individually removed and the curve reevaluated to see
% if the removal improves the ability of the gradient function to predict
% the structure of the curve based on minimising the CHI^2 error between
% prediction and observed.
% In addition the curves are analysed with a short-time fourier transform
% to visualise the frequency distribution across the timecourse of the
% experiments. Time windows associated with wider frequency distributions
% than the background curve are associated with noise. The selection of
% points for removal is weighted based on the inclusion or not of that
% point in a region of wider frequency distribution.

% Inputs    [Cexp] Unfiltered metabolite measurements across time
%           [Cexpl] List of metabolite names, natural language
%           [time] Experiment measurement timestamps
%           [pre_gradient] toggle for initialising the value of moving mean
%           based on the first [pre_gradient_range] points
%           [exp_weight] weighting factor for exponential moving average,
%           affects the relative weighting of proxmial and distal points to
%           fine tune the behaviour of the filter
%           [Exp_Err] Tolerable error between predicted and observed points
%           in moving mean calculations, ordinarily 10% but may be replaced
%           with a more involved function if desired

% Outputs   [Clean_Met] Updated Version of Cexp containing the cleaned
%           metabolite  concentrations



Fs=10;



% Importing raw data and interpolating through missing points 
% Potential to replace with the refined method seen in LoadData
for j=1:size(Cexp,2) 
    interp_met=Cexp(:,j);       
    missing=find(isnan(interp_met));            
    if isempty(missing)==0 
        if sum(missing)==length(interp_met)
            % leave array empty, will return blank in plot
        else
    t_missing=time(missing);
    interp_met(:,2)=time;
    interp_met=rmmissing(interp_met);
    interp_fn=griddedInterpolant(interp_met(:,2),interp_met(:,1));
    interp_vals=interp_fn(t_missing);
    Cexp(missing,j)=interp_vals;
        end
    end
     clear interp_met missing t_missing interp_fn interp_vals
end     





figure(10)
Clean_Met=NaN(size(Cexp,1),size(Cexp,2));

for j=1:size(Cexp,2) 
    
    
    clearvars -except j Exp_Err Clean_Met Cexp Cexpl time pre_gradient exp_weight Fs pre_gradient_range 
    
    if sum(isnan(Cexp(:,j)))~=size(Cexp,2) % if the metabolite has no recordings then it is skipped
    
    
    pt_track=0; % Tracking the number of points removed to prevent over-removal
    [Freq_curve,T_freq,lin_time]=freq_analysis(Cexp,time,Fs,j); % Performing Frequency domain analysis
    Freq_flags=round((T_freq(islocalmin(Freq_curve,'MinProminence',5))).*10); 
    flag_time(1,:)=lin_time(Freq_flags); % aligning the centers of STFT windows with linearised time increments
    flag_time(2,:)=0;
    flag_time(3,:)=5; % width of error margin, 1/2 the window width used in freq_analysis
   
 
    
    ExMet=Cexp(:,j);    % importing and constructing ExMet from Cexp
    subplot(6,7,j)
    plot(time,ExMet(:,1),'-o g') % plotting the raw data   
    title(Cexpl(j))
    hold on
   errorbar(flag_time(1,:),flag_time(2,:),flag_time(3,:),'horizontal','* r')  % plotting regions of high frequency distribution along with errorbars
        
    % in order to evaluate all points for outlier status the sequence is
    % checked in both forward and reverse directions, therefore a mirrored
    % sequence is created 
    
    
    % Constructing the ExMet matrix 
    % ExMet structure: 
    % Col 1: Metabolite Concentration , Col 2: Sample times
    % Col 3: Mirrored Concentration , Col 4: Mirrored times
    % Col 5: Linear Index , Col 6: Numeric Gradient, Col 7: Mirrored Gradient
    % Col 8: Mirrored Index
    ExMet(:,2)=time;
    ExMet(:,3)=flip(ExMet(:,1));
    ExMet(1,4)=0;
    En=size(ExMet,1);
    for t=2:size(ExMet,1)       
       ExMet(t,4)=ExMet(t-1,4)+(ExMet(En+2-t,2)-ExMet(En+1-t,2)); 
    end
    ExMet(:,5)=[1:1:size(ExMet,1)];
    ExMet(:,6)=gradient(ExMet(:,1),ExMet(:,2));    
    ExMet(:,7)=gradient(ExMet(:,3),ExMet(:,4));    
    ExMet(:,8)=flip(ExMet(:,5));
    
    % Input option: uses the mean of the first N datapoints in a sequence
    % which can account for noisy initial segments.
    if pre_gradient==1
ExMet(1,6)=mean(ExMet(1:pre_gradient_range,6));
ExMet(1,7)=mean(ExMet(1:pre_gradient_range,7));
    end
    
    % Starting scan and removal process
    Exit_flag=0;
     while Exit_flag==0
         clear BiDir_sum
        pt_track=pt_track+1;
        
        Dir=1; % scanning forwards
        [Summary,Rem_flag]=Det_Grad_Analysis(ExMet,Dir,Exp_Err,exp_weight);
        Fwd_scan=Rem_flag;
        Fwd_summary=Summary;
        clear Rem_flag
        clear Summary
        Dir=2; % scanning in reverse
        [Summary,Rem_flag]=Det_Grad_Analysis(ExMet,Dir,Exp_Err,exp_weight);
        Rev_scan=Rem_flag;
        Rev_summary=Summary;

        % Loc_weight will apply a weighting factor based on the appearance
        % of a point in the forward and/or reverse scans. Points that
        % appear in only one evaluation are negatively weighted for removal
        % Loc_list is a list of flagged locations for removal
        % candidates
        clear Loc_list Loc_weight
        Loc_weight(1,:)=ExMet(:,5);
        Loc_weight(2,:)=0;
                
        Loc_list=Fwd_scan(1,(find(Fwd_scan(2,:)==1))); Loc_list=cat(2,Loc_list,Rev_scan(1,(find(Rev_scan(2,:)==1))));
        
        for pos=1:length(Loc_list)
            P_loc=find(Loc_weight(1,:)==Loc_list(pos));
            Loc_weight(2,P_loc)=Loc_weight(2,P_loc)+1;              
        end        
        % Done constructing the Loc_weight array        
        
        if isempty(Loc_list)==1 
            Exit_flag=1;
            Rem_list=[];
            
        else % Constructing the BiDir_sum matrix which contains the
             % results of the error analysis, structure is:
             % Col 1 & 3 - no. of removal candidates flagged in fwd and rev
             % directions
             % Col 2 & 4 - Sum of CHI^2 errors in forward and reverse
             % directions
             % Col 5 - linear index of point removed

            Loc_list=(unique(Loc_list))';
            BiDir_sum(1,1:2)=Fwd_summary; BiDir_sum(1,3:4)=Rev_summary;
            BiDir_sum(1,5)=0;   % State of curve in base state

            % Outputs a list of removal candidates with associated error
            % quantification.
        for i=1:length(Loc_list)
            [Rem_sum]=Candidate_removal(ExMet,Exp_Err,i,Loc_list,exp_weight);   
            BiDir_sum(i+1,:)=Rem_sum;
        end

            % ranking removal candidates
            Baseline=BiDir_sum(1,:);
            CHI_th=Baseline(2)+Baseline(4);
            Flag_th=Baseline(1)+Baseline(3);
            BiDir_sum(1,:)=[];
            [Pt_rem]=candidate_ranking(BiDir_sum,flag_time,CHI_th,Flag_th,time,Loc_weight);
            
            
            
            Rem_list(pt_track)=Pt_rem(1,3);
            % Removing the target step from the metabolite 
            X=find(ExMet(:,5)==Pt_rem(1,3)); Y=find(ExMet(:,8)==Pt_rem(1,3));
            
            Fwd_ExMet(:,1:2)=ExMet(:,1:2); Fwd_ExMet(:,3:4)=ExMet(:,5:6);
            Rev_ExMet(:,1:2)=ExMet(:,3:4); Rev_ExMet(:,3:4)=ExMet(:,7:8);
            Fwd_ExMet(X,:)=[]; Rev_ExMet(Y,:)=[];
            clear ExMet
            ExMet(:,1:2)=Fwd_ExMet(:,1:2); ExMet(:,3:4)=Rev_ExMet(:,1:2);
            ExMet(:,5:6)=Fwd_ExMet(:,3:4); ExMet(:,7:8)=Rev_ExMet(:,3:4);
            
            if size(Pt_rem,1)==2 && pt_track < (size(Cexp,1)/5) % applied if a second removal candidate was also selected
                pt_track=pt_track+1;
                         Rem_list(pt_track)=Pt_rem(2,3);
            % Removing the target step from the metabolite 
            X=find(ExMet(:,5)==Pt_rem(2,3)); Y=find(ExMet(:,8)==Pt_rem(2,3));
            
            Fwd_ExMet(:,1:2)=ExMet(:,1:2); Fwd_ExMet(:,3:4)=ExMet(:,5:6);
            Rev_ExMet(:,1:2)=ExMet(:,3:4); Rev_ExMet(:,3:4)=ExMet(:,7:8);
            Fwd_ExMet(X,:)=[]; Rev_ExMet(Y,:)=[];
            clear ExMet
            ExMet(:,1:2)=Fwd_ExMet(:,1:2); ExMet(:,3:4)=Rev_ExMet(:,1:2);
            ExMet(:,5:6)=Fwd_ExMet(:,3:4); ExMet(:,7:8)=Rev_ExMet(:,3:4);
               
            end
            
            % exit if the removal does not significantly improve the error
            % or 25% of points have been removed
            if Pt_rem(1)==0 || (Pt_rem(2)/CHI_th)<=0.9 || pt_track>=(size(Cexp,1)/5)
                Exit_flag=1;
            end
            
            
        end
    
        
     end
     
     % exports the cleaned data and plots the cleaned set 
     Clean_Met(:,j)=Cexp(:,j); 
     Clean_Met(Rem_list,j)=NaN;
     
    else
       Clean_Met(:,j)=Cexp(:,j); 
    end
     
      plot(time,Clean_Met(:,j),'-* b')
     fprintf('Done with %i \n',j)
end

end

function [Summary,Rem_flag]=Det_Grad_Analysis(ExMet,Dir,Exp_Err,exp_weight)
% Rem_flag has the format: 
% Row 1: index location % Row 2: removal flag (1/0) % Row 3: CHIsq error
A=[1 2 5 6 ; 3 4 8 7]; % Locations of column elements in the ExMet matrix associated with forward or reverse directionality
Grad_av=0; % starting the weighted average gradient
Rem_flag(1,:)=ExMet(:,A(Dir,3))'; % indexing for clarity
for x=1:size(ExMet,1)-1

     W=exp_weight/(x+1); % updating weighting factor for the gradient function
    % Generating the new moving average gradient for this position
    Grad_av=(ExMet(x,A(Dir,4))*W)+(Grad_av*(1-W)); 
    % Making the prediction of the subsequent point based on previous
    % gradient mean
    Pt_Pred=ExMet(x,A(Dir,1))+((ExMet(x+1,A(Dir,2))-ExMet(x,A(Dir,2)))*Grad_av);
    
    % When an observed point falls outside the error window of expectation
    % the point becomes flagged for removal, in either case the CHI^2
    % distance is recorded.
    if (Pt_Pred*(1-Exp_Err))<=ExMet(x+1,A(Dir,1)) && ExMet(x+1,A(Dir,1))<=(Pt_Pred*(1+Exp_Err))
        Rem_flag(2,x+1)=0;
        Rem_flag(3,x+1)=((ExMet(x+1,A(Dir,1))-Pt_Pred)^2)/abs(Pt_Pred);
    else
        Rem_flag(2,x+1)=1;
        Rem_flag(3,x+1)=((ExMet(x+1,A(Dir,1))-Pt_Pred)^2)/abs(Pt_Pred); 
    end
    % Export the results of analysis to relevant layer
    Summary(1,1)=sum(Rem_flag(2,:));
    Summary(1,2)=sum(Rem_flag(3,:));
    
    
    
end


end


function [Rem_sum]=Candidate_removal(ExMet,Exp_Err,i,Loc_list,exp_weight)
% Removes the candidate location and re-assesses using the gradient
% analysis, outputting a summary of resulting error statistics
i=Loc_list(i);
ExMet_Base=ExMet;

% for forward scan 
ExMet(find(ExMet(:,5)==i),:)=[];
Dir=1;
[Summary,~]=Det_Grad_Analysis(ExMet,Dir,Exp_Err,exp_weight);
Fwd_Sum=Summary;

% for reverse scan
ExMet=ExMet_Base;
ExMet(find(ExMet(:,8)==i),:)=[];
Dir=2;
[Summary,~]=Det_Grad_Analysis(ExMet,Dir,Exp_Err,exp_weight);

% assemble summary of state after removal
Rev_Sum=Summary;
Rem_sum(1,1:2)=Fwd_Sum; Rem_sum(1,3:4)=Rev_Sum; Rem_sum(1,5)=i;


end

function [Pt_rem]=candidate_ranking(BiDir_sum,flag_time,CHI_th,Flag_th,time,Loc_weight)
% Ranks removal candidates according to minimisation of CHI^2 error
% Points that are flagged only in one directional assessment are negatively
% weighted (increasing error)
% Points that are located inside a window flagged as having increased
% frequency distribution are positively weighted. Fewer frequency flagged
% regions are more strongly associated with regions of noise, so when fewer
% flags are present they provide a greater weighting factor.


for i=1:size(BiDir_sum,1)
    freq_modifier=1;
    for L=1:size(flag_time,2) % applying frequency weighting
        if (flag_time(1,L)-flag_time(3,L))<=time(BiDir_sum(i,5)) && time(BiDir_sum(i,5))<=(flag_time(1,L)+flag_time(3,L))
            freq_modifier=freq_modifier-(1/size(flag_time,2));
        else
        end
    end
    
%     Loc_index=(Loc_weight(1,BiDir_sum(i,5)));
    % Loc_I_W applies a weighting factor based on the agreement of forward
    % and reverse scans
    Loc_I_W=1+(2*(1-(Loc_weight(2,find(Loc_weight(1,:)==BiDir_sum(i,5)))/2)));
    % Comparative array: Col 1 - ratio of removal candidates vs the
    % baseline, values of <1 represent a positive improvements.
    % Col 2 - weighted sum of CHI^2 errors, weighting the frequency or
    % Location modifier can adjust the results.
    Comp_array(i,1)=(BiDir_sum(i,1)+BiDir_sum(i,3))/Flag_th;
%        Comp_array(i,2)=(BiDir_sum(i,2)+BiDir_sum(i,4))*Loc_I_W;
    Comp_array(i,2)=(BiDir_sum(i,2)+BiDir_sum(i,4))*(freq_modifier)^2*Loc_I_W;
    Comp_array(i,3)=BiDir_sum(i,5);    
end

% Sort by weighted CHI^2 error
Comp_array=sortrows(Comp_array,[1 2]);

% Sends the highest ranked point for removal, if the 2nd point is close to
% the first then it will be sent for removal as well
Pt_rem(1,:)=Comp_array(1,:);

if size(Comp_array,1)>1
if Comp_array(2,2)<Comp_array(1,2)*1.25
    Pt_rem(2,:)=Comp_array(2,:);
end
end

end

function [Freq_curve,T_freq,lin_time]=freq_analysis(Cexp,time,Fs,j)

Raw_input=Cexp(:,j);
Met_start=Raw_input(1);

% linearising the input data to allow fourier analysis
for i=1:length(Raw_input)-1
   t1=time(i); t2=time(i+1); t_no=(t2-t1);
   lin_time_seg=linspace(t1,t2,t_no);
   Met1=Raw_input(i); Met2=Raw_input(i+1);
   lin_met_seg=linspace(Met1,Met2,t_no);
   
     
   if i==1
      lin_time=lin_time_seg;
      lin_met=lin_met_seg;
   else
       lin_time(end)=[];
       lin_met(end)=[];
       lin_time=cat(2,lin_time,lin_time_seg);
       lin_met=cat(2,lin_met,lin_met_seg);   
   end   
    

end

% Performing Short Time Fourier Transform, standard window size is 15 units
% this may be scaled based on the size of the linearised dataset
[s,f,t]=stft(lin_met,Fs,'Window',hann(15,'periodic'));
% STFT provides data mirrored around 0 Hz, extracting uni-directional
% results
S=s(1:8,:);
F=f(1:8,:);

% putting the frequency components into a log scale relative to the largest
% magnitude
for p=1:size(S,2)
   Freq_Mag(:,p)=(abs(S(:,p)));       
end
Top=max(Freq_Mag);
Top=max(Top);

for p=1:size(S,2)
   for r=1:size(S,1)
      Z(r,p)=20*log10(abs(S(r,p))/Top);
   end
   
   % Generating a weighted sum of frequency component distribution,
   % frequencies closer to 0 Hz are low-weighted while more distant ones
   % are high-weighted, more negative values are therefore associated with
   % higher deviation from the background frequency indicating local events
   % that may be noise/error in the data
   Freq_curve(p)=sum((Z(:,p)).*F);
end

T_freq=t;
end
