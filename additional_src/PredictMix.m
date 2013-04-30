% predict the chromatogram and binning of a uniform mixture
% seqInt - the sequence database
% seqs - the sequence indices for the mix
% returns:
% tfullSeq - the binned mixture
% chromdat - the mixed chromatogram
function [tfullSeq,chromdat]=PredictMix(seqInt,seqs)
peakwidths=0.3;
peakwidthe=0.15;
binning=1;
[tfullSeq,chromdat]=PredictFullSeq(seqInt(seqs(1),:),0,0,peakwidths,peakwidthe,binning);
tfullSeq=double(tfullSeq);
for a=2:length(seqs)
    [ttfullSeq,tchromdat]=PredictFullSeq(seqInt(seqs(a),:),0,0,peakwidths,peakwidthe,binning);
    tfullSeq=tfullSeq+double(ttfullSeq);
    chromdat=chromdat+tchromdat;
end
figure;
plot(chromdat);
title(['predicted mix']);
tfullSeq=squeeze(tfullSeq);