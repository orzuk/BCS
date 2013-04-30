% To prepare the data:
% load the sequence database
load('..\RevSeq\Reduced\Red-Rev-Name');
load('..\RevSeq\Reduced\Red-Rev-Seq');
% load the chromatogram
chromnvh=scfread('3603_mixnvh-1510R.scf');
% predict the mixture chromatogram
[predbinchromnvh,mixchrom]=PredictMixSqrVar(SeqRed2B,[13294,7852,5259,4808,4583],0.0004);
% and bin it
[binchromnvh,distdatnvh]=TestBinning(predbinchromnvh,chromnvh,1,3303,3000,150,650);

% load the predicted binning database
load('PredictedBinning.mat');