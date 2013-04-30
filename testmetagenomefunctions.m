% A script for testing various functions of the project.

sigma=0.0; % currently no noise
I = 100; % number of sequences in the mixture
%% load('R_seqIntUnique.mat'); % Warning: Heavy file to load
load('R_seqInt.mat'); % Warning: Heavy file to load
opt_method =  'matlab_linprog'; % 'greedy'
PreProcessStuff=0;


NumSpecies = size(R_seqInt,1);
SeqLen = size(R_seqInt,2);

P = randperm(NumSpecies);

if(PreProcessStuff)
    % Divide it to blocks due to memory problems
    for i=1:14
        i_is = i
        R_seq_block = R_seqInt((i-1)*10000+1:min(i*10000,NumSpecies),:);
        save( ['R_seq_block_' num2str(i)], 'R_seq_block');
    end

    % Do unique
    TotalUniqueInds = [];
    for i=1:14
        i_is = i
        load( ['R_seq_block_' num2str(i)]);
        [UniqueSeqsBlock I{i} J{i}] = unique(R_seq_block, 'rows');
        save( ['UniqueSeqsBlock' num2str(i)], 'UniqueSeqsBlock');
        TotalUniqueInds = [TotalUniqueInds I{i}'+(i-1)*10000];
    end

    DoubleInds = [];
    for i=1:size(R_seqInt,1)-2
        if(isequal(R_seqInt(i,:), R_seqInt(i+2,:)))
            DoubleInds = [DoubleInds i];
        end
        if(mod(i,1000) == 0)
            i_is = i
        end
    end

    UniqueSeqsBlock = PackUint8Seqs(R_seqInt);
    UniqueSeqsBlock10 = UniqueSeqsBlock(1:10,:);
    R_seqInt10 = UnPackUint8Seqs(UniqueSeqsBlock10);

    True_R_seqInt10 = R_seqInt(1:10,:);
    True_R_seqInt10-R_seqInt10
    figure; imagesc(    True_R_seqInt10-R_seqInt10); colorbar;


    R_seqInt = UnPackUint8Seqs(UniqueSeqsBlock);
    save 'R_seqIntUnique.mat' 'R_seqInt';

    % Do union to gain a little
    for i=1:7
        R1 = load( ['UniqueSeqsBlock' num2str(i*2-1)]);
        R2 = load( ['UniqueSeqsBlock' num2str(i*2)]);
        UniqueSeqsBlock = union(R1.UniqueSeqsBlock, R2.UniqueSeqsBlock, 'rows');
        save( ['BigUniqueSeqsBlock' num2str(i)], 'UniqueSeqsBlock');
    end


    % Count number of 'unique' seqs
    NumUniqueSeqs = 0;
    for i=1:14
        i_is = i
        load( ['UniqueSeqsBlock' num2str(i)]);
        NumUniqueSeqs = NumUniqueSeqs + size(UniqueSeqsBlock,1)
    end

    % count again
    NumUniqueSeqs = 0;
    for i=1:7
        i_is = i
        load( ['BigUniqueSeqsBlock' num2str(i)]);
        NumUniqueSeqs = NumUniqueSeqs + size(UniqueSeqsBlock,1)
    end

    % Pack 3 nucleotides in one uint8
    for i=1:7
        i_is = i
        load( ['BigUniqueSeqsBlock' num2str(i)]);
        UniqueSeqsBlock = mod(UniqueSeqsBlock, 5); % make 5 as zero - and stick to it!
        UniqueSeqsBlock = PackUint8Seqs(UniqueSeqsBlock);
        save( ['PackedBigUniqueSeqsBlock' num2str(i)], 'UniqueSeqsBlock');
    end

    % Unite everything together
    load( ['PackedBigUniqueSeqsBlock' num2str(1)]);
    BigUniqueSeqs = UniqueSeqsBlock;
    for i=2:7
        i_is = i
        load( ['PackedBigUniqueSeqsBlock' num2str(i)]);
        BigUniqueSeqs = union(BigUniqueSeqs, UniqueSeqsBlock, 'rows');
    end
    save 'PackedBigUniqueSeqs.mat' 'BigUniqueSeqs';

    % Now transfer it back to 'R_seq' format - just as unique
    NumUniqueSpecies = size(BigUniqueSeqs,1);
    SeqLen = size(BigUniqueSeqs,2)*3;

    R_seqIntUnique = PackUint8Seqs(BigUniqueSeqs);
    save( 'R_seqIntUnique.mat', 'R_seqIntUnique');

    clear R_seqInt;

    RunNum = 1; % index for running number





    DirtyHandsAndLookAtData(R_seqInt); % do all kinds of plots and look at data

else
    RunNum=1;
    % Here's the original simulation:
    [orig_ind orig_alpha Y]=SimulateMetagenomeMixture(R_seqInt,I, sigma);% generate random sequence

    if(min(min(R_seqInt)) == 0)
        R_seqInt = mod(R_seqInt+4, 5)+1; % Make sure the 'missing' is five
    end
    [rec_ind rec_alpha delta_score] = ReconstructMixture(Y, R_seqInt, I, RunNum, opt_method); % Call our reconstruction algorithm (currently greedy)
    [numofoverlap overlapset overlapalpha overlapbeta]=CompareSets(R_seqInt,orig_ind,orig_alpha,rec_ind,rec_alpha); % Compare orig&reconstructed



    %% % Get some measures of how well we did
    %% ReconstructionError = norm(alpha_vec - alpha_vec_estimate)
    %% figure; plot(alpha_vec, alpha_vec_estimate, '.'); title('True freqs. and our estimation');
    %% xlabel('True Freqs'); ylabel('Our Estimation');



end












