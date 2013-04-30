% Compare two sets 
function [numofoverlap totdist percentnotright percentmissed firstnotright firstmiss] = ...
    CompareSets6(R_seqInt,orig_ind,orig_alpha,rec_ind,rec_alpha)
% numofoverlap - the amount of sequences in both orig and reconstruction
% totdist - distance from the original data as a fraction of total original
% data
% percentnotright - fraction of incorrect reconstructed sequences
% percentmissed - fraction of missed sequences

numofoverlap=0;
overlapset=[];
unique_orig=[];
overlapalpha=[];
overlapbeta=[];

soffset=2E-6;

% the cutoff for similar sequences - 100 is to p~= 0.1%
overlap_dist=1;

% reduce the sets for homologs
%[orig_ind,orig_alpha]=ReduceSet(R_seqInt,orig_ind,orig_alpha,overlap_dist);
%[rec_ind,rec_alpha]=ReduceSet(R_seqInt,rec_ind,rec_alpha,overlap_dist);

unique_rec=[1:length(rec_ind)];

for a=1:length(orig_ind)
    orig_is_unique=1;
    for b=1:length(rec_ind)
        if SeqDist(R_seqInt,orig_ind(a),rec_ind(b))<overlap_dist
            % find sequences which do not appear in the reconstruction/orig
            % set
            orig_is_unique=0;
            rec_pos=find(unique_rec(:)==b);
            if (length(rec_pos)>0)
                unique_rec(rec_pos)=[];
            end
            numofoverlap = numofoverlap+1;
            existsinoverlap=0;
            for c=1:length(overlapset)
                if (SeqDist(R_seqInt,overlapset(c),rec_ind(b))<overlap_dist)
                    overlapalpha(c)=overlapalpha(c)+rec_alpha(b);
                    existsinoverlap=1;
                end
            end
            if (existsinoverlap==0)
                overlapset=[overlapset rec_ind(b)];
                overlapalpha=[overlapalpha rec_alpha(b)];
                overlapbeta=[overlapbeta orig_alpha(a)];
            end
        end
    end
    if (orig_is_unique)
        unique_orig=[unique_orig a];
    end
end

figure;
cla;
plot(overlapalpha,overlapbeta,'.');

tt=max(max(overlapalpha),max(overlapbeta));
xlabel('rec_alpha');
ylabel('orig_alpha');
hold on;
plot([soffset tt],[soffset tt]);
plot(zeros(length(unique_orig),1)+soffset,orig_alpha(unique_orig),'.');
unique_orig_freq=orig_alpha(unique_orig);
unique_rec_freq=rec_alpha(unique_rec);
plot(rec_alpha(unique_rec),zeros(length(unique_rec),1)+soffset,'.');
plot(sum(rec_alpha(unique_rec)),sum(orig_alpha(unique_orig)),'o');

% the total distance between observed and original frequencies
totdist=0;
% difference in correct sequences
totdist=totdist+sum((overlapalpha-overlapbeta).^2);
% add missed/incorrect sequences
totdist=totdist+sum(orig_alpha(unique_orig).^2);
totdist=totdist+sum(rec_alpha(unique_rec).^2);
totdist=sqrt(totdist);
totdist=totdist/sum(abs(orig_alpha));

percentnotright=sum(abs(rec_alpha(unique_rec)))/sum(abs(rec_alpha));
percentmissed=sum(abs(orig_alpha(unique_orig))/sum(orig_alpha));

firstnotright=max(rec_alpha(unique_rec));
if (isempty(firstnotright))
    firstnotright=0;
end
firstmiss=max(orig_alpha(unique_orig));
if (isempty(firstmiss))
    firstmiss=0;
end

disp(['totdist=' num2str(totdist)]);
title(['overlap=' num2str(numofoverlap) ' totdist=' num2str(totdist)]);
