function [k_min,E_th]=exomets_kmeans_thresholding(E_th,S,Cexp , Cexpl, time, Vars)


warning('off','SPLINES:CHCKXYWP:NaNs')

% Evaluation of steady state regions of exometabolic concentrations by
% k-mean grouping, variant for establishing the minimum number of k
% centroids
% Andrei Ligema
% 2019 DTU



Exp_err=Exomets_exp_err;

% [Cexp , Cexpl, time, Vars]=data_import;



%E_th=0.90;
for X=1:length(E_th)

for N=1:S




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
    SD_S_Err(1,j)=mean(std(S_Err));
    if ((Exp_err(j,1)+Exp_err(j,2))/(SD_S_Err(1,j)))>=E_th || (SD_S_Err(1,j))==0 
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
fprintf('Completed Iteration %u \n',X)

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