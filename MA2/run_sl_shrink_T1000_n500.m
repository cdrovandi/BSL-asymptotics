function run_sl_shrink_T1000_n500(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_1000_100reps.mat');
cov_rw = [0.00316861029133247,-0.00140450123651660;-0.00140450123651660,0.00276340540777115];


M=10000;
n=500;
K=20;
gamma = 0;


[theta, loglike] = bayes_sl_ma_acf_warton(y(:,r),M,n,cov_rw,[0.6 0.2],K,gamma);


save(['results_sl_shrink_T1000_n500_' num2str(r) '.mat'],'theta','loglike');


%delete(p);

