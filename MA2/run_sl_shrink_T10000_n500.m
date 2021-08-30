function run_sl_shrink_T10000_n500(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.000376857054321824,-0.000149134393611412;-0.000149134393611412,0.000302161744459672];

M=10000;
n=500;
K=20;
gamma = 0;


[theta, loglike] = bayes_sl_ma_acf_warton(y(:,r),M,n,cov_rw,[0.6 0.2],K,gamma);


save(['results_sl_shrink_T10000_n500_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

