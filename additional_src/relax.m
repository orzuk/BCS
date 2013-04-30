% Finaly optimization stage: Greedily Remove and add sequences to improve score.
function [rec_ind rec_alpha delta_score]=Relax(Y,R_seq,CurrentBestInds,alpha_vec,RelaxationIterations,NumToRemove,RunNum)

L = size(R_seq,2) % Sequence length
delta_score=[];
M = length(R_seq)

cY=Y;
cY(5,:)=ones(1,L);
cYY=cY;
CurError=L*4;

%relaxation
for a=1:RelaxationIterations %Relaxation iterations - remove and re-add
    % remove a creature
    olderror=CurError;
    oldInds=CurrentBestInds;
    oldalpha=alpha_vec;

    for c=1:NumToRemove
        remind=randint(1,1,length(CurrentBestInds))+1;
        rembact=CurrentBestInds(remind);
        CurrentBestInds(remind)=[];
        alpha_vec(remind)=[];
    end

    % Create and solve linear system
    X = zeros(4*L,length(CurrentBestInds));
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
    ltz=-eye(length(CurrentBestInds));
    ltzy=zeros(length(CurrentBestInds),1);
    utzy=ones(length(CurrentBestInds),1);
    oldpoint=[alpha_vec];

    %    alpha_vec = X\YY
    if (length(CurrentBestInds)==1)
        alpha_vec=min(nYY);
    else
        %        alpha_vec=linprog(-ones(cur_I,1),nX,nYY);
        options = optimset('LargeScale', 'off');
        alpha_vec=linprog(-ones(length(CurrentBestInds),1),nX,nYY,[],[],ltzy,utzy,oldpoint,options);
    end
    % Update Y vector
    cY(1:4,:) = cYY(1:4,:)- reshape(X*alpha_vec, L, 4)';

    for c=1:NumToRemove
        %    creatures_inds = CurCreatureInds(1:CurM);
        minDat=FindMaxMin(length(R_seq),size(R_seq,2),R_seq,cY);
        best_score=minDat(2);
        best_ind=minDat(1);
        %    best_ind=creatures_inds(minDat(1));

        % Add best creature
        CurrentBestInds = [CurrentBestInds  best_ind];
        % remove creature from the creature list
        %    aa=find(CurCreatureInds==best_ind);
        %    CurCreatureInds(aa)=CurCreatureInds(CurM);
        %    CurM=CurM-1;

        % Create and solve linear system
        X = zeros(4*L,length(CurrentBestInds));
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
        ltz=-eye(length(CurrentBestInds));
        ltzy=zeros(length(CurrentBestInds),1);
        utzy=ones(length(CurrentBestInds),1);

        oldpoint=[alpha_vec;0];

        %    alpha_vec = X\YY
        if (length(CurrentBestInds)==1)
            alpha_vec=min(nYY);
        else
            %        alpha_vec=linprog(-ones(cur_I,1),nX,nYY);
            options = optimset('LargeScale', 'off');
            alpha_vec=linprog(-ones(length(CurrentBestInds),1),nX,nYY,[],[],ltzy,utzy,oldpoint,options);
        end
        % Update Y vector
        cY(1:4,:) = cYY(1:4,:)- reshape(X*alpha_vec, L, 4)';
    end

    CurError = sum(sum(cY(1:4,:)));

    aa=rand(1);
    useold=(olderror<CurError);
    if (useold)
        if (aa>(a/RelaxationIterations))
            useold=0;
            useold=1;
        end
    end
    if (useold)
        [num2str(a) '- dont use - Removed ' num2str(rembact) ' from ' num2str(olderror) ' to ' num2str(CurError)]
        CurrentBestInds=oldInds;
        alpha_vec=oldalpha;
        CurError=olderror;
    else
        [num2str(a) '- use - Removed ' num2str(rembact) ' from ' num2str(olderror) ' to ' num2str(CurError)]
    end
    delta_score=[delta_score CurError];
    rec_ind=CurrentBestInds;
    rec_alpha=alpha_vec;
end
save (['Relaxation-' num2str(RunNum)],'Y','rec_ind','rec_alpha','delta_score');
