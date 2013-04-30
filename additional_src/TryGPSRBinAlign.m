% try GPSR with chromatogram binning predicted for each sequence
% gets the binned chromatogram for each sequence and the 
% binned measurement (from TestBinning)
%
% Input: 
% RSeqBin - Predicted binned chromatogram data for each sequence in the DB
% cseq - the binned experiemental chromatogram
% startpos - binned prediction start
% endpos - binned prediction end
% tau - the GPSR l1/l2 cost balance parameter
% w - coefficient enforcing sparsity for the slack variables 
% 
% Output: 
% res - result: fraction of mixture for each species. res(i) is how much of creature i is in the mixture 
% slack_vars - part of the result corresponding to 4*L slack variables 
%
function [res slack_vars] =TryGPSRBinAlign(RSeqBin,cseq,startpos,endpos,tau,w)

old_working = 0; % which method to use (add extra slack variables) 
slack_vars = []; 

% use only part of the sequence.
RSeqBin=RSeqBin(:,startpos:endpos,:);
cseq=cseq(startpos:endpos,:);

L=size(RSeqBin,2); % length of sequence 
N=size(RSeqBin,1); % number of sequences (creatures)

%% make the mix matrix
M = zeros(4*L,N); % Mixing matrix (not sparse). M(i,j) 
YY = zeros(4*L,1); % measurements YY(i) average 'A' quantity in position i 
for jj=1:4 % loop on bases 
    tt=RSeqBin(:,:,jj);  % RSeqBin(i,j,k) is creature i in position j in nucleotide k (k can be 'A', 'C', 'G', 'T', 'N')
    M(((jj-1)*L+1):(jj*L),:) = tt'; % Fill the lines of M which corrspond to nucleotide jj
    YY(((jj-1)*L+1):(jj*L)) = cseq(:,jj); % cseq(i,j) is height of nucleotide j in position i
end
if( (~old_working) && (w ~= 0) ) % New: change to enable the mix matrix incorporate errors (add a few variables)
    M = [M w*eye(4*L)]; % add a unity vector scaled by w 
end

[res]=GPSR_BB(YY,M,tau,'Debias',1,'Verbose',0); %% solve

if(~old_working)
    slack_vars = res(N+1:end); % see the additional slack variables
    ignored_eqs_inds = vec2row(find(slack_vars))
    num_ignored_eqs = length(find(slack_vars))
    res = res(1:N); % remove slack variable 
end
res=res./sum(res);

