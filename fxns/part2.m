function Bi=part2(tj,to2,i)
% if to2(i+3)-to2(i+2)==0
%     Bi=(((to2(i+4)-tj)*((tj-to2(i+1))^2))/((to2(i+4)-to2(i+1))*(to2(i+3)-to2(i+1))*(to2(i+2)-to2(i+1))))+(((tj-to2(i))/(to2(i+3)-to2(i)))*((((tj-to2(i))*(to2(i+2)-tj))/((to2(i+2)-to2(i))*(to2(i+2)-to2(i+1))))+(((to2(i+3)-tj)*(tj-to2(i+1)))/((to2(i+3)-to2(i+1))))));
% else
if to2(i+2)-to2(i+1)~=0
    Bi=(((to2(i+4)-tj)*((tj-to2(i+1))^2))/((to2(i+4)-to2(i+1))*(to2(i+3)-to2(i+1))*(to2(i+2)-to2(i+1))))+(((tj-to2(i))/(to2(i+3)-to2(i)))*((((tj-to2(i))*(to2(i+2)-tj))/((to2(i+2)-to2(i))*(to2(i+2)-to2(i+1))))+(((to2(i+3)-tj)*(tj-to2(i+1)))/((to2(i+3)-to2(i+1))*(to2(i+2)-to2(i+1))))));
else
    Bi=0;
end
    %end