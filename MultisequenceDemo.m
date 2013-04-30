% Load the reduced (distance<2 joined together), primer aligned, data
place = 'PC'; % 'WEIZMANN' or 'UNIX'
switch place
    case 'PC'
        data_dir = '\\oxygen\seq_orzuk\zuk_weizmann\research\compressed_sensing\metagenomics\data\';
    case 'UNIX'
        data_dir = '/seq/orzuk/zuk_weizmann/research/compressed_sensing/metagenomics/data/';
    case 'WEIZMANN'
        data_dir = pwd;
end
    
load(fullfile(data_dir, 'Red-Rev-Seq'));

disp('Real sequence database');

% select a subset so GPSR will work faster so only 5000 sequences
udata=randperm(size(SeqRed2B,1));
USeq=SeqRed2B(udata(1:5000),:);

% create a random mixture with 0.15 noise, 40 creatures and uniform
% distribution (flag=0)
% oi is the indices of the non-zero creatures, oa is the frequencies, y is
% the mixture
[oi oa y]=SimulateMetagenomeMixture(USeq,0.15,1,40,0);

% reconstruct the frequency vector using GPSR and only 500 nucleotides.
% tau=10.
% note that the positive requirement on the frequencies is not inserted
% (don't know how to do it with GPSR...)
res=TryGPSR(USeq(:,100:600),y(:,100:600));

% compare the original and reconstructed frequencies. See the function for
% all the values is can return.
CompareSets6(USeq(:,:),oi,oa,find(res~=0),res(res~=0));


disp('Random sequence database');
% Generate a random sequence database
RSDB=GenerateRandSeqDB(USeq);
% again create a random mixture
[oir oar yr]=SimulateMetagenomeMixture(RSDB,0.15,1,40,0);
% reconstruct it
resr=TryGPSR(RSDB(:,100:600),yr(:,100:600));
% and compare
CompareSets6(RSDB(:,:),oir,oar,find(resr~=0),resr(resr~=0));