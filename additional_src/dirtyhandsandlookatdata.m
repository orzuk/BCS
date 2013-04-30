function R=DirtyHandsAndLookAtData(R_seqInt); % do all kinds of plots and look at data


% First perform unique
[UniqueSeqs I J] = unique(R_seqInt, 'rows');

R=UniqueSeqs;


% Find the rank of the matrix
% BinarySeqs = SeqsToBinary(UniqueSeqs);

SeqsRank = rank(BinarySeqs); % See what dimension is spanned by the sequences
