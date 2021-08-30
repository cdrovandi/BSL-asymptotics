function run_sl_T10000_K5(r)
% 10000 observations and 5 autocovariances

%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.000157244082770195	6.22720196676529e-05;
6.22720196676529e-05	9.27296131650419e-05];

M=10000;
n=500;
K=5;


[theta, loglike] = bayes_sl_ma_acf(y(:,r),M,n,cov_rw,[0.6 0.2],K);


save(['results_sl_T10000_K5_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

