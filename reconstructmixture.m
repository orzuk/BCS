% Main reconstrubtion funcion for mixture chromatogram
% Try several L1 relaxation methods - e.g. linear programming
% We need to check which one is faster: L1magic package or just
% Matlab's linprog method
%
% Input:
% Y - The mixture vector
% R_seq - The sequences matrix
% RunNum - Index saying the number of the current run
%
% Output:
% rec_ind - Indexes of the reconstructed species
% rec_alpha - Frequencies of the reconstructed species
% delta_score - The improvement in the score
% opt_method - Optimization method (greedy, linprog, l1magic)
%
function [rec_ind rec_alpha delta_score] = ReconstructMixture(Y, R_seq, I, RunNum, opt_method)

M = size(R_seq,1) % Number of Species
L = size(R_seq,2) % Sequence length
CurrentBestInds = [];
delta_score=[];
alpha_vec=[];


cY=Y; cY(5,:)=ones(1,L); cYY=cY;

CurCreatureInds=[1:M];
CurM=M;

switch opt_method

    case 'greedy'
        %	do the stuff that Amnon wrote

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
                alpha_vec=linprog(-ones(cur_I,1),nX,nYY,[],[],ltzy,utzy,oldpoint,options); % first we put the inequalities

            end
            alpha_vec'
            % Update Y vector
            cY(1:4,:) = cYY(1:4,:)- reshape(X*alpha_vec, L, 4)';
            CurError = sum(sum(cY(1:4,:)))
            delta_score=[delta_score CurError];

            rec_ind=CurrentBestInds;
            rec_alpha=alpha_vec;

        end



    case 'matlab_linprog'
        SpeciesInds = 1:200;
        [Z A b A_eq b_eq] = PrepareLinearProgram(R_seq, SpeciesInds, Y);
        
        oldpoint=[alpha_vec;0];         oldpoint = zeros(1,length(Z));

        %        alpha_vec=linprog(-ones(cur_I,1),nX,nYY);
        options = optimset('LargeScale', 'off');
        
        % First one is objective function. Second ones are equality constrains.
        % Third ones are inequality constrains. Oldpoint is a guessed starting point.
        alpha_vec=linprog(Z,A,b,A_eq(1:10,:),b_eq(1:10,:),[],[],[],options);

    case 'l1magic'

    oldpoint = A'*b; % initial guess
    xp = l1eq_pd(oldpoint, A, [], b, 0.001); % call the l1magic routine

end % Switch optimization method




save ([opt_method '-' num2str(RunNum) '-I-' num2str(cur_I)],'Y','rec_ind','rec_alpha','delta_score');



% This auxillary function prepares the matrices for the linear program
% which is used to solve the L1 minimization problem
function [Z A b A_eq b_eq] = PrepareLinearProgram(R_seq, SpeciesInds, Y)

M = length(SpeciesInds); % Number of creatures participating
N = size(R_seq, 2); % Number of nucleotides
S = Uint8SeqsToBinaryMat( R_seq(SpeciesInds,:) );
YY = reshape(Y, 1, 4*N);

Z = [zeros(1,M)  -ones(1,M) ]'; % Optimization vector: 0 for x, 1 for u.

% Inequality constrains
% alpha>0 conditions
A = zeros(3*M, 2*M);
A(1:M,1:M) = -eye(M);
A(M+1:2*M,1:M) = eye(M);
A(2*M+1:3*M,1:M) = eye(M);
A(1:M,M+1:2*M) = eye(M);
A(M+1:2*M,M+1:2*M) = eye(M);

b = [-ones(1,M)./M zeros(1,2*M)]';

% Equality constrains
A_eq = [S' zeros(4*N, M)]; % Prepare room for the A's and the U's
A_eq = [A_eq' [ones(1,M) zeros(1,M)]']'; % Add the sum zero constaint
b_eq = (YY - sum(S) ./ M)';
b_eq(end+1) = 0; % the sum zero constaint

