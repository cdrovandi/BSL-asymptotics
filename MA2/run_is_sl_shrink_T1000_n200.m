function run_is_sl_shrink_T1000_n200(r)


%c = parcluster('local')
%c.JobStorageLocation = tempdir;
%p = parpool(c,16);


load('data_1000_100reps.mat');
load('results_sl_shrink_T1000_n500.mat');

M=10000;
n=200;
K=20;
gamma = 0;


[theta_is, w_is] = is_sl_ma_acf_warton(y(:,r),M,n,mean(theta(:,:,r)),2*cov(theta(:,:,r)),K,gamma);


save(['results_is_sl_shrink_T1000_n200_' num2str(r) '.mat'],'theta_is','w_is');


%delete(p);

