% Run the GPSR algorithm with chromatogram binning predicted for each sequence
% gets the binned chromatogram for each sequence and the 
% binned measurement (from BinChromatogram)
%
% Input: 
% RSeqBin - Predicted binned chromatogram data for each sequence in the DB
% cseq - the binned experiemental chromatogram
% startpos - binned prediction start
% endpos - binned prediction end
% tau - the GPSR l1/l2 cost balance parameter
%
% Output: 
% res - resulting reconstruction 
% 
function res = GPSRChromatogram(RSeqBin,cseq,startpos,endpos,tau)

% use only part of the sequence.
RSeqBin=RSeqBin(:,startpos:endpos,:);
cseq=cseq(startpos:endpos,:);

L=size(RSeqBin,2);
N=size(RSeqBin,1);

%% make the mix matrix
X = zeros(4*L,N);
YY = zeros(4*L,1);
for jj=1:4
    tt=RSeqBin(:,:,jj);
    X(((jj-1)*L+1):(jj*L),:) = tt';
    YY(((jj-1)*L+1):(jj*L)) = cseq(:,jj);
end

%% solve
[res]=GPSR_BB(YY,X,tau,'Debias',1,'Verbose',0);
% Normalize to sum 1
res=res/sum(res);
