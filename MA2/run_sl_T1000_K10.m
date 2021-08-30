function run_sl_T1000_K10(r)
% 1000 observations and 10 autocovariances

%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.00137612338695324,0.000477923334713194;0.000477923334713194,0.00112616980221712];

M=10000;
n=500;
K=10;
y = y(1:1000,r);


[theta, loglike] = bayes_sl_ma_acf(y,M,n,cov_rw,[0.6 0.2],K);


save(['results_sl_T1000_K10_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

