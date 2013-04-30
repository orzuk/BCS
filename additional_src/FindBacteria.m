% Go over all names in bactNames and find indices of names containing name
function [idx]=FindBacteria(name,bactNames)
k = strfind(bactNames, name);
a=cellfun(@isempty,k);
idx=find(a==0);
