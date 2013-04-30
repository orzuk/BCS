% Calculate the fit score between the sequence and the reference matrix
function score = calcseqscore(seq,refseq)

minscore=1;
ll=size(refseq,2);
for a=1:ll
    if (seq(a)>0)
        minscore=min(minscore,refseq(seq(a),a));
    end
end
score=minscore;
