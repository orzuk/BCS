% Run the GPSR algorithm 
function [res]=TryGPSR(RSeq,cseq)

L=size(RSeq,2);
N=size(RSeq,1);

X = zeros(4*L,N);
YY = zeros(4*L,1);
for j=1:4
    X(((j-1)*L+1):(j*L),:) = (RSeq == j)';
    YY(((j-1)*L+1):(j*L)) = cseq(j,:);
end
%[res]=GPSR_Basic(YY,X,10);
[res]=GPSR_BB(YY,X,10,'Debias',1,'Verbose',0);
res=res/sum(res);
