BCS
---
Bacterial Compressed Sensing

The package contains source and data for the Bacterial Compressed Sensing (BCS) algorithm, 
presented in the paper
"Bacterial Community Reconstruction Using Compressed Sensing"
A. Amir and O. Zuk. Journal of Computational Biology, 18(11):1723-1741 (2011) 

Please cite the above paper if using the package. 

All files are detailes in files.txt
in order to run a simulation:
1. load the greengenes database from Red-Rev-Seq.mat (can also load the names from Red-Rev-Name.m
2. create a random mixture using SimulateMetagenomeMixture.m 
3. Reconstruct it using GPSRSimulation.m 
4. check the reconstruction using SensSel2.m 

in order to reconstruct an experimental chromatogram:
1. load the predicted and binned database from PredBinSqrt.mat
or alternatively
	load a sequence database, and predict the binned chromatograms using PredictFullSeqSqrVar.m 
2. load an experimental chromatogram (such as 3603_mix5-1510R.scf)
3. Bin it using BinChromatogram.m for various offsets (if offset unknown)
4. SQRT transform both datasets (binned database and binned experimental chromatogram)
5. Reconstruct for each offset using GPSRChromatogram.m
6. Find the best offset for reconstruction (minimal reconstruction distance from original chromatogram)
Note that if the original mixture is known, offset can also be found using TestBinning for various offsets

For any questions and comments, please contact: Amnon Amir, amnonim@gmail.com or 
						Or Zuk, orzuk@ttic.edu


