function run_sl_T100_K10(r)
% 100 observations and 10 autocovariances

%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.0264809972836268,-0.00191528405859929;-0.00191528405859929,0.0319452422735003];

M=10000;
n=500;
K=10;
y = y(1:100,r);


[theta, loglike] = bayes_sl_ma_acf(y,M,n,cov_rw,[0.6 0.2],K);


save(['results_sl_T100_K10_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

