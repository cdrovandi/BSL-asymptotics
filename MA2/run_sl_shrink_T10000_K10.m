function run_sl_shrink_T10000_K10(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.000376857054321824,-0.000149134393611412;-0.000149134393611412,0.000302161744459672];

M=10000;
n=50;
K=10;
gamma = 0;
y = y(:,r);


[theta, loglike] = bayes_sl_ma_acf_warton(y,M,n,cov_rw,[0.6 0.2],K,gamma);


save(['results_sl_shrink_T10000_K10_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

