% Compute RMSE/sensitivity/selectivity statistics for various thresholds
%
% Input:
% seq - the sequence database
% oi - original mixture (non-zero frequency) species indices
% oa - original mixture (non-zero frequency) species frequencies
% res - the result of the BCS algorithm reconstruction
%
% Output:
% each variable is a vector for various thresholds
% othr - the thresholds used
% oprec - Precision
% osens - Sensitivity
% ospec - Specificity
% oprecf - frequency corrected Precision (not used in article)
% osensf - frequency corrected Sensitivity (not used in article)
% ospecf - frequency corrected Specificity (not used in article)
% otdist - RMSE
%
function [othr,oprec,osens,ospec,oprecf,osensf,ospecf,otdist]= ...
    SensSel2(seq,oi,oa,res)
thr=1E-4;
othr=[];
oprec=[];
osens=[];
ospec=[];
oprecf=[];
osensf=[];
ospecf=[];
otdist=[];

olen=length(oi);
%thefig=figure;
mr=max(res);
while (thr<=0.4)
    disp(thr);
    cres=find(res>thr);
    [numofoverlap totdist percentnotright percentmissed firstnotright firstmiss]=CompareSets5(seq,oi,oa,cres,res(cres),[]);
    othr=[othr; thr];
    osens=[osens; numofoverlap/olen];
    ospec=[ospec; (size(seq,1)-length(cres)-length(oi)+numofoverlap)/(size(seq,1)-length(oi))];
    oprec=[oprec; numofoverlap/length(cres)];
    otdist=[otdist; totdist];
    
    osensf=[osensf; 1-percentmissed];
    oprecf=[oprecf; 1-percentnotright];
    ospecf=[ospecf; (size(seq,1)-length(cres)-length(oi)+numofoverlap)/(size(seq,1)-length(oi))];
    thr=thr*1.5;
end
