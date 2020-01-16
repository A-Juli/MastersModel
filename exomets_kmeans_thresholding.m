function [k_min,E_th]=exomets_kmeans_thresholding(E_th,S,Cexp_clean, Cexpl, time, Vars_clean)


warning('off','SPLINES:CHCKXYWP:NaNs')
% Andrei Ligema
% 2019 DTU/CFB
% Evaluation of steady state regions of exometabolic concentrations by
% k-mean grouping, variant for establishing the minimum number of k
% centroids
% Assembles data from across replicates to ensure that k is common for each
% set of experiments
%
% Inputs    [E_th] sensitivity for convergence between experimental data
%           and reconstruction based on cluster algorithm. Reccomended at
%           90% (0.9) but may be lowered to reduce number of clusters or
%           raised to increase.
%           [S] Number of repeats, k-means is a stochastic process so the
%           confidence for number of clusters will improve with more
%           [Cexp_clean] Table of experimental measurements
%           [Cexpl] Metabolite names in natural language
%           [time] timestamps of measurements
%           [Vars_clean] experiment metadata
% see exomets_process and exomets_data_import for additional information
% regarding inputs and the procedure
%
% Outputs   [k_min] minimum number of clusters required to achieve error
%           convergence between data and reconstruction
%           [E_th] unused
Cexpl_base=Cexpl;
X=1;    

% evaluating variation in experimental observations for all metabolites
% standard deviation is taken per timepoint per metabolite across
% replicates
% the mean of this deviation per metabolite is then used as the threshold
% for convergence

for p=1:size(Cexpl,2) 
    clear Meas_err Meas
    for s=1:size(Cexp_clean,2)
    Meas(:,s)=Cexp_clean{1,s}{1,1}(:,p);    
    end
    Meas=rmmissing(Meas);
   for r=1:size(Meas,1)
       Meas_err(r)=std(Meas(r,:));
       %Meas_err(r)=abs((((Meas(r)*1.1)-Meas(r))^2)/Meas(r));      
   end
   Meas_err=rmmissing(Meas_err);
   if isempty(Meas_err)==1
      Meas_err=0; 
   end
   Exp_err(p,1)=mean(Meas_err);
   Exp_err(p,2)=std(Meas_err);
end
Exp_err=cat(1,Exp_err,Exp_err);


% k-means is based on q for metabolites, need to construct the appropriate
% inputs here, in this procedure all replicates are run together in order
% to provide a single k value


q=[];
Interp_Cexp=[];
for T=1:size(Cexp_clean,2)
    
 Cexp=Cexp_clean{1,T}{1,1};
 Vars=Vars_clean{1,T}{1,1};
    
[q_rep,Met_Exp_rep,Delta_Vars_rep,Interp_Cexp_rep,Interp_Vars,xq]=k_means_input(time,Cexp,Vars);

q=cat(2,q,q_rep);
Interp_Cexp=cat(2,Interp_Cexp,Interp_Cexp_rep);

end

if size(Cexp_clean,2)>=2
for T=2:size(Cexp_clean,2)
   Cexpl=cat(2,Cexpl,Cexpl_base);
end
end


% start of procedure proper, runs for S repeats to narrow 95% CI of k


for N=1:S






k=1;
CHI_Flag=0;

while CHI_Flag==0

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

clear errcount
for j=1:size(Ret_C,2)
    clear S_Err
    S_Err(:,1)=Interp_Cexp(:,j);
    S_Err(:,2)=Ret_C(:,j);
    S_Err=S_Err';
    Fit_Err=abs(((S_Err(2,:)-S_Err(1,:)).^2)./S_Err(1,:));
    Fit_Err=rmmissing(Fit_Err);
    SD_S_Err(1,j)=mean(Fit_Err);
    if isnan(SD_S_Err(1,j))==1
    SD_S_Err(1,j)=0;
    end
    if (((Exp_err(j,1)+Exp_err(j,2)))/(SD_S_Err(1,j)))>=E_th || (SD_S_Err(1,j))==0 
    errcount(j)=1;
    else
        errcount(j)=0;
    end
end
   


for j=1:size(Ret_C,2)
    Err_prop(k,j)=SD_S_Err(1,j)/Exp_err(j,1);
end


sum(errcount);
clear Err
% end of error convergence process

if k>=2
%    used if evaluating error threshold
if  sum(errcount)==size(Ret_C,2)%abs(CHI_diffs(k-1))<=abs((CHI_diffs_max))
    CHI_Flag=1;
    
else
    k=k+1;
end
else
    k=k+1;
end


end

K_rec(X,N)=k;



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



clear output_q



end




CI95=tinv([0.975],((size(K_rec,2))-1));
for i=1:length(E_th)
    K_av(i)=mean(K_rec(i,:));
    K_SEM(i)=std(K_rec(i,:))/sqrt(size(K_rec,2));
    K_CI95(i)=K_SEM(i)*CI95;
end
figure(3)
errorbar(E_th,K_av,K_CI95)
xlabel('Error Threshold');
ylabel('mean no. of k')

for i=1:length(E_th)
    k_min(i)=ceil(K_av(i)+K_CI95(i));
end

end

function [q_rep,Met_Exp_rep,Delta_Vars_rep,Interp_Cexp_rep,Interp_Vars,xq]=k_means_input(time,Cexp,Vars)

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
        Delta_Vars_rep(i,j)=Interp_Vars(i+1,j)-Interp_Vars(i,j);
    end
end
% produces a table of delta values for time, VCD, vol and biomass
% done establishing bioreactor variables


for j=1:size(Cexp,2)
    Met_Exp_rep(:,1)=time;
splinefun=csaps(time,Cexp(:,j),p);

for i=1:length(xq)
    Interp_Cexp_rep(i,j)=ppval(splinefun,xq(i));
end

    

end
% generates a full table of interpolated concentration values for each
% metabolite



for j=1:size(Interp_Cexp_rep,2)
    for i=1:length(xq)-1
        DeltaC(i,j)=Interp_Cexp_rep(i+1,j)-Interp_Cexp_rep(i,j);
    end
end
% full table of differences for the metabolites


for j=1:size(DeltaC,2)
    for i=1:size(DeltaC,1)
        q_rep(i,j)=DeltaC(i,j)/Interp_Vars(i,4);
    end
end
% evaluate q for the table of differences

% done preparing the inputs for kmeans

end
