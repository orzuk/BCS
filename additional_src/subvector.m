function [aa,maxv]=subvector(a,b)
% find the best occurance of a in b (and the match length)
ll=length(a);
score=zeros(ll,1);
for aa=1: length(b)+1-ll
    score(aa)=sum(b(aa:aa+ll-1)==a);
end
[maxv aa]=max(score);
