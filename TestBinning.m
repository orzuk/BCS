% Compare two binned chromatograms
% MixBin - the binned predicted chromatogram
% seqbin - the binned experimental chromatogram
% returns:
% difdat - the sqr difference
% difdist - the difference values
function [difdat,difdist]=TestBinning(MixBin,seqbin,startpos,endpos)
seqbin=seqbin/mean(sum(seqbin(startpos:endpos,:),2));
MixBin=MixBin/mean(sum(MixBin(startpos:endpos,:),2));

minlen=min(length(seqbin(:,1)),length(MixBin(:,1)));
difdat=(seqbin(1:minlen,1)-MixBin(1:minlen,1)).^2;
difdat=difdat+(seqbin(1:minlen,2)-MixBin(1:minlen,2)).^2;
difdat=difdat+(seqbin(1:minlen,3)-MixBin(1:minlen,3)).^2;
difdat=difdat+(seqbin(1:minlen,4)-MixBin(1:minlen,4)).^2;

difdist=seqbin(startpos:endpos,1)-MixBin(startpos:endpos,1);
difdist=[difdist;seqbin(startpos:endpos,2)-MixBin(startpos:endpos,2)];
difdist=[difdist;seqbin(startpos:endpos,3)-MixBin(startpos:endpos,3)];
difdist=[difdist;seqbin(startpos:endpos,4)-MixBin(startpos:endpos,4)];
