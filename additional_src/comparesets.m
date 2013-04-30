% Compare two sets&scores - evaluate performance of the algorithm's output
% by comparing it to the true mixture.
function [numofoverlap overlapset overlapalpha overlapbeta]=CompareSets(R_seqInt,orig_ind,orig_alpha,rec_ind,rec_alpha)
% randres=[];
% randpoints=randint(5000,2,[1,length(R_seqInt)]);
% for a=1:5000
%     randres=[randres SeqDist(R_seqInt,randpoints(a,1),randpoints(a,2))];
% end
% distdist=randres;
% hist(randres,500);

numofoverlap=0;
overlapset=[];
unique_orig=[];
overlapalpha=[];
overlapbeta=[];

% the cutoff for similar sequences - 100 is to p~= 0.1%
overlap_dist=20;

% reduce the sets for homologs
[orig_ind,orig_alpha]=ReduceSet(R_seqInt,orig_ind,orig_alpha,overlap_dist);
[rec_ind,rec_alpha]=ReduceSet(R_seqInt,rec_ind,rec_alpha,overlap_dist);

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
plot(overlapalpha,overlapbeta,'.');
tt=max(max(overlapalpha),max(overlapbeta));
xlabel('rec_alpha');
ylabel('orig_alpha');
hold on;
plot([0 tt],[0 tt]);
plot(zeros(length(unique_orig),1),orig_alpha(unique_orig),'x');
plot(rec_alpha(unique_rec),zeros(length(unique_rec),1),'+');
plot(sum(rec_alpha(unique_rec)),sum(orig_alpha(unique_orig)),'o');

