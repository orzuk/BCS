% Compare reconstructed mixture to original mixture and calculate statistics
%
% Input:
% R_seqInt - the sequence database
% orig_ind - indices of species present in the original mixture
% orig_alpha - Their frequencies
% rec_ind - indices of the species present in the reconstruction
% rec_alpha - Their frequencies
% thefig - figure handle where to plot the results
%
% Output:
% numofoverlap - the amount of sequences in both orig and reconstruction
% totdist - distance from the original data as a fraction of total original data
% percentnotright - fraction of incorrect reconstructed sequences
% percentmissed - fraction of missed sequences
% firstnotright - the frequency of the largest incorrect species
% firstmiss - the frequency of the highest missed species
%
function [numofoverlap totdist percentnotright percentmissed firstnotright firstmiss] = ...
    CompareSets5(R_seqInt,orig_ind,orig_alpha,rec_ind,rec_alpha,thefig)
numofoverlap=0;
overlapset=[];
unique_orig=[];
overlapalpha=[];
overlapbeta=[];

overlap_dist=1;

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

unique_orig_freq=orig_alpha(unique_orig);
unique_rec_freq=rec_alpha(unique_rec);
if (~isempty(thefig))
    figure(thefig);
    cla;
    plot(overlapalpha,overlapbeta,'.');

    tt=max(max(overlapalpha),max(overlapbeta));
    xlabel('rec_alpha');
    ylabel('orig_alpha');
    hold on;
    plot([0 tt],[0 tt]);
    plot(zeros(length(unique_orig),1),orig_alpha(unique_orig),'x');
    plot(rec_alpha(unique_rec),zeros(length(unique_rec),1),'+');
    plot(sum(rec_alpha(unique_rec)),sum(orig_alpha(unique_orig)),'o');
end

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
