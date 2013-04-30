% find the closest bacteria in the db to the sequence
% seqred2b - the sequence database used for the reconstruction (non-binned)
% seq - the sequence to compare
function [res]=FindDBSeq(SeqRed2B,seq,startpos,endpos)
res=[];
for a=1:size(SeqRed2B,1)
res=[res length(find(SeqRed2B(a,startpos:endpos)~=seq(startpos:endpos)))];
end
[retval,ret]=min(res);
disp(retval);
disp(ret);

