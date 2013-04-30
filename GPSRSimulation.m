% Perform reconstruction of a simualted mixture using GPSR
%
% Input:
% RSeq - The sequence database (N creatures, L nucleotides) (a=1,c=2,g=3,t=4)
% cseq - the PSSM of the simulated mixture (L nucleotides * 4)
%
% Output:
% res - the (sparse) frequency vector of each of the database species
%
function res = GPSRSimulation(RSeq,cseq)

L=size(RSeq,2);
N=size(RSeq,1);

X = zeros(4*L,N);
YY = zeros(4*L,1);
for j=1:4
    X(((j-1)*L+1):(j*L),:) = (RSeq == j)';
    YY(((j-1)*L+1):(j*L)) = cseq(j,:);
end

tau=10;
[res]=GPSR_BB(YY,X,tau,'Debias',1,'Verbose',0);
% normalize to sum 1
res=res/sum(res);
