% Here is the reconstrubtion funcion. Current greedy algorithm sucks.
% We assume that we know I. If not, we can enumerate on it and give some
% penalty in order to choose the optimal I
function [rec_ind rec_alpha delta_score] = ReconstructMixture(Y, R_seq, I, RunNum)


L = size(R_seq,2) % Sequence length
CurrentBestInds = [];
delta_score=[];
alpha_vec=[];

M = length(R_seq)

cY=Y;
cY(5,:)=ones(1,L);
cYY=cY;

CurCreatureInds=[1:M];
CurM=M;

for cur_I=1:I %loop on creatures - each time add one

    % All creatures
%    creatures_inds = CurCreatureInds(1:CurM);
    minDat=FindMaxMin(length(R_seq),size(R_seq,2), R_seq, cY);
    best_score=minDat(2);
    best_ind=minDat(1);
%    best_ind=creatures_inds(minDat(1));
    
    % Add best creature
    CurrentBestInds = [CurrentBestInds  best_ind]
    % remove creature from the creature list
%    aa=find(CurCreatureInds==best_ind);
%    CurCreatureInds(aa)=CurCreatureInds(CurM);
%    CurM=CurM-1;
    
    % Create and solve linear system
    X = zeros(4*L,cur_I);
    YY = zeros(4*L,1); 
    for j=1:4
        X(((j-1)*L+1):(j*L),:) = (R_seq(CurrentBestInds,:) == j)';
        YY(((j-1)*L+1):(j*L)) = cYY(j,:);
    end    
    % remove equations which mean nothing (due to early sequence termination)
    trueEQ=find(sum(X,2)>0);
    nX=X(trueEQ,:);
    nYY=YY(trueEQ);
    
    % add the alpha>0 conditions
    ltz=-eye(cur_I);
    ltzy=zeros(cur_I,1);
    utzy=ones(cur_I,1);
%    nX=[nX ; ltz];
%    nYY=[nYY ; ltzy];
    oldpoint=[alpha_vec;0];

%    alpha_vec = X\YY
    if (cur_I==1)
        alpha_vec=min(nYY);
    else
%        alpha_vec=linprog(-ones(cur_I,1),nX,nYY);
        options = optimset('LargeScale', 'off');
        alpha_vec=linprog(-ones(cur_I,1),nX,nYY,[],[],ltzy,utzy,oldpoint,options);
    end
    alpha_vec'
% Update Y vector
    cY(1:4,:) = cYY(1:4,:)- reshape(X*alpha_vec, L, 4)';
    CurError = sum(sum(cY(1:4,:)))
    delta_score=[delta_score CurError];
    
    rec_ind=CurrentBestInds;
    rec_alpha=alpha_vec;
    
end

save (['Greedy-' num2str(RunNum) '-I-' num2str(cur_I)],'Y','rec_ind','rec_alpha','delta_score');
