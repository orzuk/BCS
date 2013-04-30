% test the performance of the reconstruction alg. for noisy binned data
% predbinchrom - the predicted chromatogram
% noise - the standard deviation of the predicted noise
function [res]=TestNoise(predbinchrom,predbinchromall,startpos,endpos,tau,noise)
noisepredbinchrom=predbinchrom+randn(size(predbinchrom))*noise;
res=TryGPSRBinAlign(predbinchromall,noisepredbinchrom,startpos,endpos,tau);