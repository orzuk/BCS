disp('Simulated MixAll');
disp('tau=16000, region=150:650');
pallT16kR150650=TryGPSRBinAlign(predbinSqr,predbinchromall,150,650,16000);
disp('tau=10000, region=150:650');
pallT10kR150650=TryGPSRBinAlign(predbinSqr,predbinchromall,150,650,10000);
disp('tau=25000, region=150:650');
pallT25kR150650=TryGPSRBinAlign(predbinSqr,predbinchromall,150,650,25000);

disp('tau=16000, region=150:350');
pallT16kR150350=TryGPSRBinAlign(predbinSqr,predbinchromall,150,350,16000);
disp('tau=25000, region=150:350');
pallT125kR150350=TryGPSRBinAlign(predbinSqr,predbinchromall,150,350,25000);

disp('Experimental MixAll');
disp('tau=16000, region=150:650');
allT16kR150650=TryGPSRBinAlign(predbinSqr,binchromall,150,650,16000);
disp('tau=10000, region=150:650');
allT10kR150650=TryGPSRBinAlign(predbinSqr,binchromall,150,650,10000);
disp('tau=25000, region=150:650');
allT25kR150650=TryGPSRBinAlign(predbinSqr,binchromall,150,650,25000);

disp('tau=16000, region=150:350');
allT16kR150350=TryGPSRBinAlign(predbinSqr,binchromall,150,350,16000);
disp('tau=10000, region=150:350');
allT10kR150350=TryGPSRBinAlign(predbinSqr,binchromall,150,350,10000);
disp('tau=25000, region=150:350');
allT25kR150350=TryGPSRBinAlign(predbinSqr,binchromall,150,350,25000);



disp('Simulated NVH');
disp('tau=10000, region=150:650');
pnvhT10kR150650=TryGPSRBinAlign(predbinSqr,predbinchromnvh,150,650,10000);
disp('tau=25000, region=150:650');
pnvhT25kR150650=TryGPSRBinAlign(predbinSqr,predbinchromnvh,150,650,25000);
disp('tau=40000, region=150:650');
pnvhT40kR150650=TryGPSRBinAlign(predbinSqr,predbinchromnvh,150,650,40000);

disp('tau=10000, region=150:350');
pnvhT10kR150350=TryGPSRBinAlign(predbinSqr,predbinchromnvh,150,350,10000);
disp('tau=25000, region=150:350');
pnvhT25kR150350=TryGPSRBinAlign(predbinSqr,predbinchromnvh,150,350,25000);
disp('tau=40000, region=150:350');
pnvhT40kR150350=TryGPSRBinAlign(predbinSqr,predbinchromnvh,150,350,40000);

disp('Experimental Mixnvh');
disp('tau=10000, region=150:650');
nvhT10kR150650=TryGPSRBinAlign(predbinSqr,binchromnvh,150,650,10000);
disp('tau=25000, region=150:650');
nvhT25kR150650=TryGPSRBinAlign(predbinSqr,binchromnvh,150,650,25000);
disp('tau=40000, region=150:650');
nvhT40kR150650=TryGPSRBinAlign(predbinSqr,binchromnvh,150,650,40000);

disp('tau=10000, region=150:350');
nvhT10kR150350=TryGPSRBinAlign(predbinSqr,binchromnvh,150,350,10000);
disp('tau=16000, region=150:350');
nvhT16kR150350=TryGPSRBinAlign(predbinSqr,binchromnvh,150,350,16000);
disp('tau=16000, region=150:350');
nvhT25kR150350=TryGPSRBinAlign(predbinSqr,binchromnvh,150,350,25000);
disp('tau=16000, region=150:350');
nvhT40kR150350=TryGPSRBinAlign(predbinSqr,binchromnvh,150,350,40000);

