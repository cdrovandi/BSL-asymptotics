function run_sl_shrink_T1000_K20(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_10000_100reps.mat');
cov_rw = [0.00371242178803702,-0.00182877177794663;-0.00182877177794663,0.00345338033284252];

M=10000;
n=50;
K=20;
gamma = 0;
y = y(1:1000,r);


[theta, loglike] = bayes_sl_ma_acf_warton(y,M,n,cov_rw,[0.6 0.2],K,gamma);


save(['results_sl_shrink_T1000_K20_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

