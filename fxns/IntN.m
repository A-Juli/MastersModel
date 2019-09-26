%author Veronica Martinez, The University of Queensland, AIBN
%19-Jan-2015

% to estimate matrix Me to fit metabolic conentrations:
% c(t)=c_0+S*K*P*Me*IN(t)
function Me=IntN(to,order)

Me=zeros(length(to)-order,length(to)-1);
for i=1:size(Me,2)
    for j=1:min(i,size(Me,1))
        Me(j,i)=(to(j+order)-to(j));
    end
end
Me=Me/order;

function Me=IntN3(to)
to=[to to(end) to(end) to(end)];
Me=zeros(length(to)-6,length(to)-4);
for i=1:size(Me,2)
    for j=1:min(i,size(Me,1))
        Me(j,i)=(to(j+3)-to(j));
    end
end
Me=Me/3;
