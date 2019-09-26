function [Q_Pred,idx,Ret_C]=exomets_kmeans(k,Cexp , Cexpl, time, Vars, spreadsheet,page)


warning('off','SPLINES:CHCKXYWP:NaNs')

% Evaluation of steady state regions of exometabolic concentrations by
% k-mean grouping, uses a threshold determined by the sister function
% exomets_kmeans_thresholding to choose the optimum number of centroids
% Andrei Ligema
% 2019 DTU
% outputs:
% idx is the cluster assignment for each time span
% C is the q value for each metabolite for each of the clusters
% sumd is sum of distances to each cluster, within cluster
% D is the distances to each centroid from each point

% retreiving the metabolite concentrations from the kmeans data
% preliminary implementation, depends on all groupings being continuous

% q values are assigned to each interval


% Exp_err=Exomets_exp_err;

% [Cexp , Cexpl, time, Vars , spreadsheet,page]=data_import;

% E_th=0.90; 
reps=5;
% parameters for weak form spline fitting, to be replaced by Veronica's
% splines
xq=linspace(time(1),time(end),500);

p=5e-4;

Interp_Vars(:,1)=xq;

for j=2:size(Vars,2)
    splinefun=csaps(time,Vars(:,j),p);
    for i=1:length(xq)
        Interp_Vars(i,j)=ppval(splinefun,xq(i));
    end
end

V=1;
for j=1:length(xq)
    if xq(j)>time(V)
        V=V+1;
    else
    end
    Interp_Vars(j,3)=Vars(V,3);
end
% replaces sample volume as a stepped function consistent with experimental
% measurements

% Interp_Vars(:,4)=Interp_Vars(:,2).*Interp_Vars(:,3).*0.36;


% Gives values of time, VCD, vol, OUR and biomass for the linspace xq

for j=1:size(Interp_Vars,2)
    for i=1:length(xq)-1
        Delta_Vars(i,j)=Interp_Vars(i+1,j)-Interp_Vars(i,j);
    end
end
% produces a table of delta values for time, VCD, vol and biomass
% done establishing bioreactor variables


for j=1:size(Cexp,2)
    Met_Exp(:,1)=time;
splinefun=csaps(time,Cexp(:,j),p);

for i=1:length(xq)
    Interp_Cexp(i,j)=ppval(splinefun,xq(i));
end

    

end
% generates a full table of interpolated concentration values for each
% metabolite



for j=1:size(Interp_Cexp,2)
    for i=1:length(xq)-1
        DeltaC(i,j)=Interp_Cexp(i+1,j)-Interp_Cexp(i,j);
    end
end
% full table of differences for the metabolites


for j=1:size(DeltaC,2)
    for i=1:size(DeltaC,1)
        q(i,j)=DeltaC(i,j)/Interp_Vars(i,4);
    end
end
% evaluate q for the table of differences

% done preparing the inputs for kmeans


% k=1;
% CHI_Flag=0;
% Q_Pred=zeros(1:length(Cexpl));
% while CHI_Flag==0

[idx,C,~,~]=kmeans(q,k);

for j=1:size(C,2)
    for i=1:size(idx,1)
Q_Pred(i,j)=C(idx(i),j);
    end
end
for j=1:size(q,2)
    for i=1:size(q,1)
        obs=abs(q(i,j));
        exp=abs(Q_Pred(i,j));
        Diffs(i,j)=((obs-exp)^2)/exp;
    end
end



% retrieval of concentrations for error measurements
% 
% Ret_DeltaC=zeros(size(idx),size(C,2));
% for j=1:size(C,2)
%     for i=1:size(idx,1)
% Ret_DeltaC(i,j)=Q_Pred(i,j)*Interp_Vars(i,4);
%     end
% end
% 
% Ret_C(1,:)=Interp_Cexp(1,:);
% for j=1:size(C,2)
%     for i=2:size(xq,2)
%         Ret_C(i,j)=Ret_C(i-1,j)+Ret_DeltaC(i-1,j);
%     end
% end
% 
% clear errcount
% for j=1:size(Ret_C,2)
%     clear S_Err
%     S_Err(:,1)=Interp_Cexp(:,j);
%     S_Err(:,2)=Ret_C(:,j);
%     S_Err=S_Err';
%     SD_S_Err(1,j)=mean(std(S_Err));
%     if ((Exp_err(j,1)+Exp_err(j,2))/(SD_S_Err(1,j)))>=E_th || (SD_S_Err(1,j))==0 
%     errcount(j)=1;
%     else
%         errcount(j)=0;
%     end
% end
%    
% 
% 
% for j=1:size(Ret_C,2)
%     Err_prop(k,j)=SD_S_Err(1,j)/Exp_err(j,1);
% end
% 
% 
% sum(errcount);
% clear Err
% end of error convergence process

% if k>=2
%     % used if evaluating error threshold
% % if  sum(errcount)==size(Ret_C,2)%abs(CHI_diffs(k-1))<=abs((CHI_diffs_max))
% %     CHI_Flag=1;
% %     
% % else
% %     k=k+1;
% % end
% else
%     k=k+1;
% end
% 
% 
% end



% Repeat best fit to adjust for stochastic variability in k-means
% choose the best version of the fit
repeats={};

for r=1:reps
[idx_r,C_r]=kmeans(q,k);
repeats{r,1}={idx_r};
repeats{r,2}={C_r};

for j=1:size(C_r,2)
    for i=1:size(idx_r,1)
        
Q_Pred_r(i,j)=repeats{r,2}{1}(repeats{r,1}{1}(i),j);
    end
end
for j=1:size(q,2)
    for i=1:size(q,1)
        obs=abs(q(i,j));
        exp=abs(Q_Pred_r(i,j));
        Diffs(i,j)=((obs-exp)^2)/exp;
    end
end
CHI_r(r)=nansum(Diffs,'all');
end
CHI_best=min(CHI_r);
CHI_best_loc=find(CHI_r==CHI_best,1);

idx=cell2mat(repeats{CHI_best_loc,1});
C=cell2mat(repeats{CHI_best_loc,2});

for j=1:size(C,2)
    for i=1:size(idx,1)
Q_Pred(i,j)=C(idx(i),j);
    end
end
for j=1:size(q,2)
    for i=1:size(q,1)
        obs=abs(q(i,j));
        exp=abs(Q_Pred(i,j));
        Diffs(i,j)=((obs-exp)^2)/exp;
    end
end

CHI2=nansum(Diffs,'all');

% Repeats the last level 5 times and returns the best fitting version
%%%%

% Retrieves concentration curves from the q predictions
for j=1:size(C,2)
    for i=1:size(idx,1)
Ret_DeltaC(i,j)=Q_Pred(i,j)*Interp_Vars(i,4);
    end
end

Ret_C(1,:)=Interp_Cexp(1,:);
for j=1:size(C,2)
    for i=2:size(xq,2)
        Ret_C(i,j)=Ret_C(i-1,j)+Ret_DeltaC(i-1,j);
    end
end

% linearises the assignment of centroid numbers for clarity
Lin=1;
IDX_lin(1)=1;
for j=2:length(idx)
    if idx(j)~=idx(j-1)
        Lin=Lin+1;
    IDX_lin(j)=Lin;
    else
        IDX_lin(j)=Lin;
    end
end


% 
fprintf('Fitting complete on sample %s %s after iteration %u \n',spreadsheet,page,k);
fprintf('Fit %u was selected from 5 repeats of the target number of knots \n',CHI_best_loc);

output_q=input('Would you like to view the fitted curves? (Y/N) \n','s');
if output_q=='Y'
    figure(1)
% plot(xq(1:499),idx','o')
% hold on
plot(xq(1:499),IDX_lin,'*')
ylabel('Group Assignment');
xlabel('Time Interval');
ylim([0,k]);
 
figure(2)
for i=1:size(Cexp,2)
    subplot(6,6,i)
    plot(Err_prop(:,i))
    title(Cexpl(i))
end

figure(3)
for i=1:size(Cexp,2)
    subplot(6,6,i)
    plot(xq,Interp_Cexp(:,i),'r')
    hold on
    plot(xq,Ret_C(:,i),'b')
    title(Cexpl(i))
end
figure(5)
subplot(1,2,1)
plot(CHI_rec)
xlabel('k')
ylabel('p')
subplot(1,2,2)
plot(CHI_diffs)
xlabel('dk')
ylabel('dp')

figure(6)
for i=1:size(Cexp,2)
    subplot(6,6,i)
    plot(Q_Pred(:,i),'r')
    title(Cexpl(i))
end




else
end


clear output_q


end

    

% outputs:
% idx is the cluster assignment for each time span
% C is the q value for each metabolite for each of the clusters
% sumd is sum of distances to each cluster, within cluster
% D is the distances to each centroid from each point

% retreiving the metabolite concentrations from the kmeans data
% preliminary implementation, depends on all groupings being continuous

% q values are assigned to each interval


% fully retrieved set of concentrations



% validation ouputs

% figure(1)
% plot(xq(1:499),idx','o')
% ylabel('Group Assignment');
% xlabel('Time Interval');
% ylim([0,k]);

% XQ=linspace(1,500,500);
% figure(2)
% subplot(1,2,1)
% plot(xq,Interp_Cexp(:,2),'r')
% hold on
% plot(xq,Ret_C(:,2),'b')
% legend('b-spline','k-mean derived')
% xlabel('time: hrs')
% ylabel('concentration')
% title('spline vs k-mean derived concentrations')
% % ylim([0,30])
% subplot(1,2,2)
% plot(time,Cexp(:,2),'-*')
% xlabel('time: hrs')
% ylabel('concentration')
% title('raw data for concentration')
% ylim([0,30])

% figure(4)
% for i=1:size(Cexp,2)
%     subplot(6,6,i)
%     plot(xq,Interp_Cexp(:,i),'r')
%     hold on
%     plot(xq,Ret_C(:,i),'b')
%     title(Cexpl(i))
% end

% figure(3)
% CT=C';
% bar(CT)
% xlabel('Metabolite no')
% ylabel('Predicted q from k-mean')
% legend('k1','k2','k3','k4')


% fprintf('Completed Iteration %u \n',X)
% 
% 
% 
% CI95=tinv([0.975],((size(K_rec,2))-1));
% for i=1:length(E_th)
%     K_av(i)=mean(K_rec(i,:));
%     K_SEM(i)=std(K_rec(i,:))/sqrt(size(K_rec,2));
%     K_CI95(i)=K_SEM(i)*CI95;
% end
% figure(3)
% errorbar(E_th,K_av,K_CI95)
% xlabel('Error Threshold');
% ylabel('mean no. of k')
