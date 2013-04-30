% What the hell is this? not used currently
function simseq=find_similar(R_seqInt,oind)
simseq=[];
orig_dat=R_seqInt(oind,:);
for (a=1:length(R_seqInt))
    if ((R_seqInt(a,:)-orig_dat)==0)
        simseq=[simseq a];
    end
end
length(simseq)
