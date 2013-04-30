% Read a fasta file 
function [seq]=LoadFasta(filename)
fseq=fastaread(filename);
ffseq=fseq.Sequence;
seq=zeros(size(ffseq));
for a=1:length(ffseq)
    if (ffseq(a)=='N')
        seq(a)=5;
    end
    if (ffseq(a)=='A')
        seq(a)=1;
    end
    if (ffseq(a)=='C')
        seq(a)=2;
    end
    if (ffseq(a)=='G')
        seq(a)=3;
    end
    if (ffseq(a)=='T')
        seq(a)=4;
    end
end
