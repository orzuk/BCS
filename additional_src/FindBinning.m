function FindBinning(seqInt,peakWidth)
[ofs]=PredictFullSeq(seqInt,1,1,peakWidth,1);
res=[];
resb=[];
for a=0.9:0.01:1.1
    fs=ofs;
    [ts]=PredictFullSeq(seqInt,1,1,peakWidth,a);
    ts=ts(size(fs));
    fs=fs(size(ts));
    b=ts-fs;
    b=b.^2;
    c=sum(sum(b));
    resb=[resb a];
    res=[res c];
end
figure;
plot(resb,res);
