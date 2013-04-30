% Bins the chromatogram (from scfread) into constant sized bins following
% height noramlization
% inputs:
% cseq - the scf structure for the experimentally measured chromatogram
% bin size - the function multiplies it by the default bin size which is
% 12, (we use 1)
% offset - The position difference between the measured and predicted chromatograms
% startpos - binned prediction start
% endpos - binned prediction end

% returns:
% seqbin - the binned measured chromatogram
function [seqbin]=BinChromatogram(cseq,binsize,offset,startpos,endpos)

%% Bin the cseq chromatogram
% the size of each bin (native)
cres=12;
nb=binsize*cres;
cfs=1;
chrmodat=zeros(length(cseq.A),4);
chromdat(:,1)=cseq.A;
chromdat(:,2)=cseq.C;
chromdat(:,3)=cseq.G;
chromdat(:,4)=cseq.T;
chromdatall=sum(chromdat,2);
chromdatall=runmean(chromdatall,250);
chromdat=chromdat./repmat(chromdatall,1,4);
for (cpos=1:nb:size(chromdat,1)-nb)
    if (round(cpos+offset)<=0)
        fullSeq(cfs,:)=[0 0 0 0];
    else
        if (round(cpos+offset+nb)>size(chromdat,1))
            fullSeq(cfs,:)=[0 0 0 0];
        else
            fullSeq(cfs,:)=sum(chromdat(round(cpos+offset):round(cpos+offset+nb),:));
        end
    end
    cfs=cfs+1;
end
seqbin=fullSeq;

seqbin=seqbin/mean(sum(seqbin(startpos:endpos,:),2));
