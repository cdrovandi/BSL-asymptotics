function run_sl_T1000_n500(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_1000_100reps.mat');
cov_rw = [0.00129844336604171,0.000663163192888927;0.000663163192888927,0.000995020153158127];

M=10000;
n=500;
K=20;


[theta, loglike] = bayes_sl_ma_acf(y(:,r),M,n,cov_rw,[0.6 0.2],K);


save(['results_sl_T1000_n500_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

