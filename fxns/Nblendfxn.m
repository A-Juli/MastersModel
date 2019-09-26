%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

% B-spline blending function of generic order
function [N, IN]=Nblendfxn(tj,to,order)

to2= [ones(1,order-1)*min(to) to ones(1,order-1)*max(to)];
total=length(to2);

if tj == max(to2)
    position = find(tj==to2,1)-1;
else
    position = find(tj>=to2(1:end-1) & tj<to2(2:end));
end

N = zeros(total-1,1);
N(position) = 1;

k = 2;
while k <=order
    N_pre = N;
    N = zeros(total-k,1);
    
    rowNZ = find(N_pre~=0)';
    if isempty(rowNZ)
        N = zeros(total-order-2,1);
        return
    end
    rowNZ = [rowNZ(1)-1 rowNZ];

    for position = rowNZ
        if to2(position+k-1)-to2(position) == 0
            A = 0;
        else
            A = (tj-to2(position))/(to2(position+k-1)-to2(position))*N_pre(position);
        end
        if to2(position+k)-to2(position+1) == 0
            B = 0;
        else
            B = (to2(position+k)-tj)/(to2(position+k)-to2(position+1))*N_pre(position+1);
        end
        N(position) =  A+B;
    end    
    k = k+1;
end


to2 = [min(to2) to2 max(to2)];
total = length(to2);
if nargout >1
    N_pre = [0; N; 0];
    IN = zeros(total-k,1);
    
    rowNZ = find(N_pre~=0)'; 
    if isempty(rowNZ)
        IN = zeros(total-order-2,1);
        return
    end
    rowNZ = [rowNZ(1)-1 rowNZ];
    
    for position = rowNZ
        if to2(position+k-1)-to2(position) == 0
            A = 0;
        else
            A = (tj-to2(position))/(to2(position+k-1)-to2(position))*N_pre(position);
        end
        if to2(position+k)-to2(position+1) == 0
            B= 0;
        else
            B= (to2(position+k)-tj)/(to2(position+k)-to2(position+1))*N_pre(position+1);
        end
        IN(position) =  A+B;
    end
end

N=N(order:end-order+1);
IN=IN(order+1:end);
