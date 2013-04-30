% Find homologus sequences and cluster them. 
%
% Input:
% 
% R_seqInt - A vector of sequences
% oind - indexes of sequences we want to check
% oalpha - alpha weights for sequences we want to check
% overlap_dist - the maximal distance below which two sequences are
% considered identical. This is any distance given by the function 'SeqDist'
%
% Output:
% rind - indexes of output sequences 
% ralpha - alpha weights for output sequences 
%
function [rind,ralpha] = ...
    ReduceSet(R_seqInt,oind,oalpha,overlap_dist)
rind=[];
ralpha=[];
for a=1:length(oind)
    found=0;
    for b=1:length(rind)
        if (SeqDist(R_seqInt,oind(a),rind(b))<overlap_dist)
            ralpha(b)=ralpha(b)+oalpha(a);
            found=1;
            break;
        end
    end
    if (~found)
        rind=[rind oind(a)];
        ralpha=[ralpha oalpha(a)];
    end
end
['reduced from ' num2str(length(oind)) ' to ' num2str(length(rind))]
