% Predict full sequence chromatogram, and divide to bins using quadratic peak position correction 
% (i.e. not linear position dependence but also a sqrVal*pos^2 correction)
% calculate total A,C,G,T in each bin
%
% Input:
% posfile - the table of the local-sequence dependent peak heights (e.g.
% Positions_Stats.txt)
% heightfile - the table of the local-sequence dependent peak heights (e.g.
% Heights_Stats.txt)
% seqInt - the dna sequences for which to draw the chromatogram
% peak width - the relative width of each peak (we use 0.4)
% binning - the bin sizes for the intensity integral (we use 1)
% sqrVal is the square streching factor value (we use sqrVal=0.00036)
%
% Output:
% fullSeq - a ~seqlen(depends on sequence and binning) on 5 (ACGT+unknown)
% predicted intensity matrix
%
function [tfullSeq,chromdat]=PredictFullSeqSqrVar(posfile,heightfile,seqInt,peakwidth,binning,sqrVal)

sLen=size(seqInt,2);

disp(['Loading Position File ' posfile]);
posDat=load(posfile);
numOfNucs=log2(size(posDat,1))/2;

disp(['Loading Height File ' heightfile]);
heightDat=load(heightfile);
if (numOfNucs~=log2(size(heightDat,1))/2)
    display('different sizes of data files!!!');
    return
end
disp(['Sequences contains ' num2str(numOfNucs) ' nucleotides (inc. end)']);

seqRepStep=100;

% init the out sequence to ones (so if not enough data - define as 1)
posVar=ones(size(seqInt,2),'single');
heightVar=ones(size(seqInt,2),'single');

% the n-mers values of each sequence position
tmpseqs=zeros(numOfNucs,sLen);
% to scan for Ns in the sequence
tseqok=zeros(numOfNucs,sLen);

tfullSeq=zeros(size(seqInt,1),round(sLen/binning)+2,5,'uint8');

%figure;
%hold on;

for cseq=1:size(seqInt,1)
    if (mod(cseq,seqRepStep)==0)
        disp(cseq);
    end

    for cpos=1:numOfNucs
        tmpseqs(cpos,numOfNucs:end)= single((seqInt(cseq,cpos:end+cpos-numOfNucs)-1)) * (4^(cpos-1));
        tseqok(cpos,numOfNucs:end)= (seqInt(cseq,cpos:end+cpos-numOfNucs)>0).*(seqInt(cseq,cpos:end+cpos-numOfNucs)<5);
    end
    % the code for each position
    totseqs=sum(tmpseqs);

    % is this position ok (all nucleotides known)
    totseqok=logical(prod(tseqok));

    % if ok fill height from table
    posVar(totseqok)=posDat(totseqs(totseqok)+1,5);
    heightVar(totseqok)=heightDat(totseqs(totseqok)+1,3);

    % else put 1
    posVar(not(totseqok))=1;
    heightVar(not(totseqok))=1;


    fullSeq=zeros(round(sLen/binning)+2,5);

    seqPos=cumsum(posVar);

    %figure;
    %plot(seqPos,heightVar);

    %% draw the chromatogram
    cres=12;
    pwid=round(cres*4*peakwidth);
    chromdat=zeros(cres*sLen+1000,5);
    for cnuc=1:sLen
        y = gaussmf([-pwid:pwid],[peakwidth*cres 0])*heightVar(cnuc);
        y=y';
        
        % the peak center position
        cenpos=(seqPos(cnuc)*cres+pwid+1);
        % the quadratic peak position correction
        cenpos=cenpos+cnuc*(cnuc-1)*sqrVal/2;
        
        cenpos=round(cenpos);
        
        chromdat([cenpos-pwid:cenpos+pwid],seqInt(cseq,cnuc)) = chromdat([cenpos-pwid:cenpos+pwid],seqInt(cseq,cnuc))+y;
    end
%    figure;plot(chromdat);
    %% now fill the binned data
    nb=binning*cres;
    cfs=1;
    for (cpos=pwid+1+seqPos(1)*cres-nb/2:nb:size(chromdat,1)-nb)
        disp(cpos);
        fullSeq(cfs,:)=sum(chromdat(round(cpos):round(cpos+nb),:));
        cfs=cfs+1;
    end
    tfullSeq(cseq,:,:)=uint8(round(fullSeq(1:round(sLen/binning)+2,:)*10));
end