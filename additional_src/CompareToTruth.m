% check the homology of an unknown seq to the database
% returns the closest species (minrespos) and its hamming distance (minres)
% from the set compset
function [minspec,minrespos,minres]=CompareToTruth(who,SeqRed2B,compset,startpos,endpos)
res=[];
for a=1:length(compset)
    res=[res length(find(SeqRed2B(compset(a),startpos:endpos)~=SeqRed2B(who,startpos:endpos)))];
%    disp([num2str(compset(a)) ' gives ' num2str(res(a))]);
end
[minres,minrespos]=min(res);
minspec=compset(minrespos);

%disp(find(SeqRed2B(compset(minrespos),50:600)~=SeqRed2B(who,50:600)));