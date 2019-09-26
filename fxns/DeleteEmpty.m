%author Veronica Martinez, The University of Queensland, AIBN
%25-Sept-2013

% to delete empty spaces
function withoutempty=DeleteEmpty(withempty)
for i = length(withempty):-1:1
    if isempty(withempty{i})
        withempty(i)=[];
    end
end
withoutempty=withempty;
    